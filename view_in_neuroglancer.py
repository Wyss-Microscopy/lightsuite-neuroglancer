# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
#     "neuroglancer==2.41.*",
#     "tifffile>=2023.7.10",
#     "imagecodecs>=2023.1.23",
#     "numpy>=1.24,<3",
#     "pandas>=2.0,<4",
#     "scipy>=1.10,<2",
#     "scikit-image>=0.21,<1",
#     "pyyaml>=6,<7",
# ]
# ///
# Keep in sync with [project].dependencies in pyproject.toml. Exact tested versions
# are in view_in_neuroglancer.py.lock (uv lock --script view_in_neuroglancer.py).
# The Python cap is deliberate: on a newer Python some packages (imagecodecs,
# scikit-image) have no prebuilt wheels yet, so pip/uv try to compile them, which
# fails on Windows without Visual C++. Raise it once they publish wheels.
"""View LightSuite registration results in Neuroglancer.

Usage:
    uv run view_in_neuroglancer.py [RESULTS_FOLDER] [options]

If RESULTS_FOLDER is omitted a folder picker opens (or you are asked in the
terminal). The folder (or any folder below it) must contain LightSuite's
volume_registered/sample_space/sample_space_manifest.json. The registration
results are converted once to Neuroglancer "precomputed" data (cached in
<results>/neuroglancer/), served locally, and opened in the browser with:

  - ch1, ch2, ...   : sample channels on the registration grid
  - template        : atlas template warped onto the sample (hidden)
  - boundaries      : atlas region outlines
  - one segmentation layer per level of the atlas hierarchy (e.g. Allen
    substructure / structure / division / category) with region names,
    colours, tags (search e.g. "#Isocortex") and per-region statistics
  - extra label volumes (e.g. spinal cord segments)
  - region meshes   : 3D surfaces of every region at every hierarchy level
  - a side panel    : search a region by name ("locus coeruleus") to highlight it and
                      jump to it, and choose a LUT (Grays, Fire, Viridis, ...) per channel
  - with --zoom     : raw multires ROI channels at native resolution, placed with
                      the saved ROI-to-overview transforms (no re-registration)

Region names / hierarchy come from, in order of preference:
  1. Fiederling spinal cord atlas Atlas_Regions.csv (parent tree)
  2. Allen parcellation_to_parcellation_term_membership.csv
  3. a Gubra / Perens LSFM ARA2_annotation_info.csv region table
  4. a downloaded BrainGlobe atlas' structures.json (~/.brainglobe,
     /home/*/.brainglobe or $BRAINGLOBE_DIR), e.g. Perens LSFM, Princeton
  5. the sample's stats/*region_stats.csv (region / structure / division)
  6. bare label numbers
Atlas files are looked for in --atlas-dir, the config's atlas.atlas_dir (and
the subfolders of both) and ~/atlases/*.
"""

import argparse
import difflib
import functools
import hashlib
import http.server
import json
import os
import re
import secrets
import shutil
import struct
import sys
import threading
import unicodedata
import urllib.parse
import webbrowser
from pathlib import Path

import neuroglancer
import numpy as np
import pandas as pd
import tifffile
import yaml
from scipy import ndimage
from skimage.measure import marching_cubes

__version__ = "0.2.0"
TOOL_VERSION = 4
MANIFEST = Path("volume_registered") / "sample_space" / "sample_space_manifest.json"
CHUNK = 64
SYNTHETIC_ID_BASE = 3_000_000_000
CHANNEL_COLORS = ["#ffffff", "#00ff00", "#ff00ff", "#00ffff", "#ffff00"]
ALLEN_MEMBERSHIP = "parcellation_to_parcellation_term_membership.csv"
CORD_REGIONS = "Atlas_Regions.csv"
REGION_TABLE = "ARA2_annotation_info.csv"  # Gubra / Perens LSFM atlas region table


# --------------------------------------------------------------------------
# Locating results
# --------------------------------------------------------------------------

def find_results_dirs(root: Path, max_depth: int = 4) -> list[Path]:
    """Folders under root (inclusive) holding a LightSuite sample_space manifest."""
    root = root.resolve()
    if root.name == "sample_space" and (root / MANIFEST.name).is_file():
        return [root.parent.parent]
    found = []
    for dirpath, dirnames, _ in os.walk(root):
        d = Path(dirpath)
        if (d / MANIFEST).is_file():
            found.append(d)
        depth = len(d.relative_to(root).parts)
        dirnames[:] = [] if depth >= max_depth else sorted(
            n for n in dirnames if not n.startswith(".") and n != "neuroglancer")
    return found


def ask_folder() -> Path:
    try:
        import tkinter
        from tkinter import filedialog
        tk = tkinter.Tk()
        tk.withdraw()
        chosen = filedialog.askdirectory(title="Select LightSuite registration results folder")
        tk.destroy()
        if not chosen:
            sys.exit("No folder selected.")
        return Path(chosen)
    except Exception:  # no display / no tkinter
        return Path(input("Path to LightSuite registration results folder: ").strip().strip("'\""))


def choose_results_dir(arg: str | None) -> Path:
    root = Path(arg).expanduser() if arg else ask_folder()
    if not root.is_dir():
        sys.exit(f"Not a folder: {root}")
    cands = find_results_dirs(root)
    if not cands:
        sys.exit(f"No LightSuite results found in {root}\n(looked for {MANIFEST})")
    if len(cands) == 1:
        return cands[0]
    print("Several registration results found:")
    for i, c in enumerate(cands, 1):
        print(f"  [{i}] {c}")
    k = input(f"Choose 1-{len(cands)}: ").strip()
    try:
        return cands[int(k) - 1]
    except (ValueError, IndexError):
        sys.exit("Invalid choice.")


# --------------------------------------------------------------------------
# Sample metadata
# --------------------------------------------------------------------------

class Sample:
    def __init__(self, root: Path):
        self.root = root
        self.ss_dir = root / MANIFEST.parent
        self.manifest = json.loads((root / MANIFEST).read_text(encoding="utf-8"))
        self.perm = self.manifest.get("permute_sample_to_atlas") or [1, 2, 3]
        self.atlas_permuted = bool(self.manifest.get("atlas_volumes_permuted", False))
        self.config = self._read_config()
        self.stats = self._read_stats("region_stats")  # measured in atlas space
        self.stats_sample = self._read_stats("region_stats_sample")  # in sample space
        self.atlas_name = self._atlas_name()
        self.name = self._sample_name()

    def _read_config(self) -> dict:
        for p in sorted(self.root.glob("*.y*ml")):
            try:
                cfg = yaml.safe_load(p.read_text(encoding="utf-8"))
            except Exception:
                continue
            if isinstance(cfg, dict) and isinstance(cfg.get("atlas"), dict):
                return cfg
        return {}

    def _read_stats(self, stem: str) -> pd.DataFrame | None:
        """stats/<stem>.csv (all channels), or else the stats/chanNN_<stem>.csv files."""
        p = self.root / "stats" / f"{stem}.csv"
        if p.is_file():
            return pd.read_csv(p)
        parts = sorted((self.root / "stats").glob(f"chan[0-9][0-9]_{stem}.csv"))
        return pd.concat([pd.read_csv(f) for f in parts]) if parts else None

    def _atlas_name(self) -> str:
        name = (self.config.get("atlas") or {}).get("provider")
        if not name and self.stats is not None and "atlas" in self.stats:
            name = self.stats["atlas"].iloc[0]
        if not name and (self.root / "regopts.json").is_file():
            name = json.loads((self.root / "regopts.json").read_text(encoding="utf-8")).get("brain_atlas")
        return str(name or "unknown")

    def _sample_name(self) -> str:
        name = (self.config.get("sample") or {}).get("name")
        if not name and self.stats is not None and "sample" in self.stats:
            name = self.stats["sample"].iloc[0]
        return str(name or self.root.name)

    def resolve(self, p: str | None) -> Path | None:
        """A manifest path, or the same file name inside the results folder if moved."""
        if not p:
            return None
        p = Path(p)
        if p.is_file():
            return p
        for d in (self.ss_dir, self.root):
            if (d / p.name).is_file():
                return d / p.name
        return None

    def atlas_dirs(self, override: str | None) -> list[Path]:
        dirs = []
        cfg_dir = (self.config.get("atlas") or {}).get("atlas_dir")
        for d in (override, cfg_dir):
            if not d:
                continue
            # the folder itself, then its subfolders: a folder holding several
            # atlases (e.g. D:\atlas\{allen_CCF3_10um, perens_20um}) also works
            d = Path(d).expanduser()
            dirs.append(d)
            try:
                dirs += sorted(s for s in d.iterdir() if s.is_dir())
            except OSError:
                pass
        if self.manifest.get("segments_csv"):
            dirs.append(Path(self.manifest["segments_csv"]).parent)
        atlases = Path.home() / "atlases"
        if atlases.is_dir():
            dirs += sorted(d for d in atlases.iterdir() if d.is_dir())
        return [d for d in dict.fromkeys(dirs) if d.is_dir()]

    def source_files(self) -> list[Path]:
        files = [self.root / MANIFEST]
        files += sorted(self.root.glob("*.y*ml"))  # config: atlas_dir, brainglobe_name
        files += [self.resolve(v) for k, v in self.manifest.items() if k.endswith("_path")]
        files += [self.resolve(v) for v in self.manifest.get("channel_paths", {}).values()]
        files += [self.root / "registration_analysis" / "region_metrics.csv",
                  self.root / "stats" / "region_stats.csv",
                  self.root / "stats" / "region_stats_sample.csv"]
        files += sorted((self.root / "stats").glob("chan[0-9][0-9]_region_stats*.csv"))
        return [f for f in files if f is not None and f.is_file()]


# --------------------------------------------------------------------------
# Volumes
# --------------------------------------------------------------------------

def load_yxz(path: Path) -> np.ndarray:
    """Multi-page TIFF as (Y, X, Z), like lightsuite's load_registration_volume."""
    with tifffile.TiffFile(path) as tf:
        return np.moveaxis(np.stack([p.asarray() for p in tf.pages]), 0, -1)


def permute_brain_volume(volume: np.ndarray, permvec: list[int]) -> np.ndarray:
    """Permute and flip axes (lightsuite permute_brain_volume / permuteBrainVolume.m)."""
    out = np.transpose(volume, [abs(v) - 1 for v in permvec])
    for dim, val in enumerate(permvec):
        if val < 0:
            out = np.flip(out, axis=dim)
    return out


def to_xyz(vol_yxz: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(vol_yxz.transpose(1, 0, 2))


def compute_boundaries(annot: np.ndarray) -> np.ndarray:
    edge = np.zeros(annot.shape, dtype=np.uint8)
    for ax in range(3):
        a = np.swapaxes(annot, 0, ax)
        e = np.swapaxes(edge, 0, ax)
        diff = a[1:] != a[:-1]
        e[1:] |= diff
        e[:-1] |= diff
    edge[annot == 0] = 0
    return edge


def load_volumes(sample: Sample):
    """All volumes as (x, y, z) arrays on the registration grid used by the pipeline.

    Channels always get permute_sample_to_atlas (as in
    load_permuted_registration_volume); atlas volumes only when the manifest
    says they were not written already permuted.
    """
    m = sample.manifest

    def atlas_vol(path):
        v = load_yxz(path)
        return v if sample.atlas_permuted else permute_brain_volume(v, sample.perm)

    channels = {}
    for ch, p in sorted(m.get("channel_paths", {}).items(), key=lambda kv: int(kv[0])):
        path = sample.resolve(p)
        if path is None:
            print(f"  warning: channel {ch} file not found: {p}")
            continue
        channels[f"ch{ch}"] = to_xyz(permute_brain_volume(load_yxz(path), sample.perm))

    annot_path = sample.resolve(m.get("annotation_path"))
    if annot_path is None:
        sys.exit("Annotation volume (annotation_path in the manifest) not found.")
    annot = to_xyz(atlas_vol(annot_path))

    images = dict(channels)
    if (p := sample.resolve(m.get("template_path"))) is not None:
        images["template"] = to_xyz(atlas_vol(p))
    if (p := sample.resolve(m.get("boundary_path"))) is not None:
        images["boundaries"] = to_xyz(atlas_vol(p))
    else:
        images["boundaries"] = compute_boundaries(annot)

    extra = {}
    for key, p in m.items():
        if not key.endswith("_path") or key in ("annotation_path", "template_path", "boundary_path"):
            continue
        path = sample.resolve(p)
        if path is None:
            continue
        v = to_xyz(atlas_vol(path))
        if v.any():
            extra[key[: -len("_path")]] = v

    for name, v in {**images, **extra}.items():
        if v.shape != annot.shape:
            sys.exit(f"Shape mismatch: {name} {v.shape} vs annotation {annot.shape}")

    # voxel_um describes the registration grid (shape_yxz) the volumes are on now
    vox = m.get("voxel_um", [20.0, 20.0, 20.0])
    voxel_nm = [float(vox[1]) * 1000, float(vox[0]) * 1000, float(vox[2]) * 1000]
    return channels, images, annot, extra, voxel_nm


# --------------------------------------------------------------------------
# Ontology: label index -> one term per hierarchy level
# --------------------------------------------------------------------------

class Ontology:
    """levels: names, finest first. chains: label -> term key per level.
    terms: key -> dict(acronym, name, color, native_id)."""

    def __init__(self, source, levels, chains, terms, level_tags=None):
        self.source, self.levels, self.chains, self.terms = source, levels, chains, terms
        self.level_tags = level_tags or [tag_name(lv) for lv in levels]
        self.dir = None  # atlas folder the names came from

    def assign_ids(self):
        ids, used = {}, set()
        for key, t in self.terms.items():
            nid = t.get("native_id")
            if nid is not None and 0 < nid < SYNTHETIC_ID_BASE and nid not in used:
                ids[key] = nid
                used.add(nid)
        nxt = SYNTHETIC_ID_BASE
        for key in self.terms:
            if key not in ids:
                ids[key] = nxt
                nxt += 1
        self.ids = ids


def unknown_term(label):
    return ("label", int(label)), {"acronym": str(label), "name": f"label {label} (not in atlas)",
                                   "color": None, "native_id": int(label)}


def complete(onto: Ontology, labels) -> Ontology:
    """Labels the provider does not know get their own term at the finest level
    and no parent (None) at the coarser ones."""
    for lab in labels:
        if lab not in onto.chains:
            key, term = unknown_term(lab)
            onto.terms[key] = term
            onto.chains[lab] = [key] + [None] * (len(onto.levels) - 1)
    return onto


def allen_ontology(csv: Path, labels) -> Ontology | None:
    memb = pd.read_csv(csv)
    known = set(memb["parcellation_index"])
    if np.mean([lab in known for lab in labels]) < 0.9:
        return None
    order = (memb.groupby("parcellation_term_set_name")["term_set_order"].first()
             .sort_values(ascending=False))
    levels = list(order.index)
    terms = {}
    for _, r in memb.drop_duplicates("parcellation_term_label").iterrows():
        m = re.fullmatch(r"AllenCCF-Ontology-\d+-(\d+)", r["parcellation_term_label"])
        terms[r["parcellation_term_label"]] = {
            "acronym": str(r["parcellation_term_acronym"]),
            "name": str(r["parcellation_term_name"]),
            "color": str(r["color_hex_triplet"]),
            "native_id": int(m.group(1)) if m else None,
        }
    chains = {}
    sub = memb[memb["parcellation_index"].isin(labels)]
    for idx, g in sub.groupby("parcellation_index"):
        by_level = dict(zip(g["parcellation_term_set_name"], g["parcellation_term_label"]))
        chain = [by_level.get(lv) for lv in levels]
        for i in range(len(chain)):  # fill gaps with the next finer term
            if chain[i] is None:
                chain[i] = chain[i - 1] if i else None
        if chain[0] is not None:
            chains[int(idx)] = chain
    return Ontology(f"Allen membership ({csv})", levels, chains, terms)


def table_records(csv: Path) -> list[dict]:
    """Region table with id, parent id, name, acronym, red, green, blue columns:
    the Fiederling spinal cord Atlas_Regions.csv (parent_ID) or the Gubra / Perens
    ARA2_annotation_info.csv (parent_id). Black means no colour."""
    df = pd.read_csv(csv, encoding="utf-8-sig")
    parent_col = next(c for c in df.columns if c.lower() == "parent_id")
    recs = []
    for _, r in df.iterrows():
        rgb = (int(r["red"]), int(r["green"]), int(r["blue"]))
        recs.append({"id": int(r["id"]), "acronym": str(r["acronym"]),
                     "name": str(r["name"]).strip(),
                     "color": "#%02x%02x%02x" % rgb if any(rgb) else None,
                     "parent": int(r[parent_col]) if pd.notna(r[parent_col]) else 0})
    return recs


def brainglobe_records(structures_json: Path) -> list[dict]:
    recs = []
    for x in json.loads(structures_json.read_text(encoding="utf-8")):
        path = x.get("structure_id_path") or [x["id"]]
        recs.append({"id": int(x["id"]), "acronym": str(x["acronym"]), "name": str(x["name"]),
                     "color": "#%02x%02x%02x" % tuple(x["rgb_triplet"]) if x.get("rgb_triplet") else None,
                     "parent": int(path[-2]) if len(path) > 1 else 0})
    return recs


def brainglobe_dirs(sample: Sample) -> list[Path]:
    """Downloaded BrainGlobe atlases, best name match first."""
    roots = [Path(os.environ["BRAINGLOBE_DIR"])] if os.environ.get("BRAINGLOBE_DIR") else []
    roots += [Path.home() / ".brainglobe"]
    # other users' atlases: /home/* (Linux), C:\Users\* (Windows), /Users/* (macOS)
    for users in dict.fromkeys([Path.home().parent, Path("/home")]):
        try:
            roots += sorted(users.glob("*/.brainglobe"))
        except OSError:
            pass
    dirs = []
    for r in dict.fromkeys(roots):
        try:
            dirs += [d for d in sorted(r.iterdir()) if (d / "structures.json").is_file()]
        except OSError:
            continue
    atlas = sample.config.get("atlas") or {}
    want = str(atlas.get("brainglobe_name") or "")
    res = atlas.get("resolution_um")
    def score(d):
        n = d.name.lower()
        return (bool(want) and n.startswith(want.lower()),
                sample.atlas_name.lower() in n,
                res is not None and f"_{int(float(res))}um" in n)
    return sorted(dirs, key=score, reverse=True)


def tree_ontology(records: list[dict], labels, source: str, min_coverage=0.0) -> Ontology | None:
    terms = {r["id"]: {"acronym": r["acronym"], "name": r["name"], "color": r["color"],
                       "native_id": r["id"]} for r in records}
    parent = {r["id"]: r["parent"] for r in records}
    if not labels or np.mean([lab in terms for lab in labels]) < min_coverage:
        return None
    paths = {}  # label -> chain from the root down to the label
    for lab in labels:
        if lab not in terms:
            continue
        chain, cur = [lab], lab
        while parent.get(cur, 0) in terms and parent[cur] not in chain:
            cur = parent[cur]
            chain.append(cur)
        paths[lab] = chain[::-1]
    # levels counted from the root so each layer is one consistent depth;
    # the finest layer is always the labels themselves
    max_len = max((len(c) for c in paths.values()), default=1)
    depths = list(range(min(max_len - 2, 3), -1, -1))
    chains = {lab: [c[-1]] + [c[min(d, len(c) - 1)] for d in depths] for lab, c in paths.items()}
    levels, level_tags = ["regions"], ["regions"]
    for k, d in enumerate(depths, 1):
        acrs = sorted({terms[c[k]]["acronym"] for c in chains.values()})
        levels.append(f"depth {d}: " + ", ".join(acrs[:3]) + (", ..." if len(acrs) > 3 else ""))
        level_tags.append(f"depth{d}")
    # parents in this atlas are black; colour them by the mean of their leaves
    for key, t in terms.items():
        if t["color"] is None:
            leaf_cols = [terms[c[0]]["color"] for c in chains.values()
                         if key in c and terms[c[0]]["color"]]
            if leaf_cols:
                rgb = np.mean([[int(h[i:i + 2], 16) for i in (1, 3, 5)] for h in leaf_cols], 0)
                t["color"] = "#%02x%02x%02x" % tuple(int(v) for v in rgb)
    return Ontology(source, levels, chains, terms, level_tags)


def stats_ontology(stats: pd.DataFrame, labels) -> Ontology:
    df = stats.drop_duplicates("parcellation_index").set_index("parcellation_index")
    name_of = dict(zip(df["acronym"], df["name"]))
    levels = ["regions", "structure", "division"]
    terms, chains = {}, {}
    for lab in labels:
        if lab not in df.index:
            continue
        r = df.loc[lab]
        keys = [("regions", int(lab))]
        terms[keys[0]] = {"acronym": str(r["acronym"]), "name": str(r["name"]),
                          "color": None, "native_id": int(lab)}
        for lv in levels[1:]:
            acr = str(r[lv]) if lv in r and pd.notna(r[lv]) else "unassigned"
            keys.append((lv, acr))
            terms[keys[-1]] = {"acronym": acr, "name": str(name_of.get(acr, acr)),
                               "color": None, "native_id": None}
        chains[int(lab)] = keys
    return Ontology("sample stats/region_stats.csv", levels, chains, terms)


def build_ontology(sample: Sample, labels, atlas_dir_override) -> Ontology:
    labels = [int(x) for x in labels if x != 0]
    dirs = sample.atlas_dirs(atlas_dir_override)
    onto = None
    if "segments_csv" in sample.manifest or sample.atlas_name == "fiederling":
        for d in dirs:
            if (d / CORD_REGIONS).is_file():
                onto = tree_ontology(table_records(d / CORD_REGIONS), labels,
                                     f"atlas region tree ({d / CORD_REGIONS})")
                onto.dir = d
                break
    prefer = list(dict.fromkeys(
        [d for d in dirs if sample.atlas_name.lower() in d.name.lower()] + dirs))
    if onto is None:
        for d in prefer:
            if (d / ALLEN_MEMBERSHIP).is_file() and sample.atlas_name in ("allen", "unknown"):
                onto = allen_ontology(d / ALLEN_MEMBERSHIP, labels)
                if onto:
                    onto.dir = d
                    break
    if onto is None:
        for d in prefer:
            if (d / REGION_TABLE).is_file():
                onto = tree_ontology(table_records(d / REGION_TABLE), labels,
                                     f"atlas region table ({d / REGION_TABLE})", min_coverage=0.9)
                if onto:
                    onto.dir = d
                    break
    if onto is None:
        for d in brainglobe_dirs(sample):
            onto = tree_ontology(brainglobe_records(d / "structures.json"), labels,
                                 f"BrainGlobe atlas ({d})", min_coverage=0.9)
            if onto:
                onto.dir = d
                break
    stats = sample.stats if sample.stats is not None else sample.stats_sample
    if onto is None and stats is not None:
        onto = stats_ontology(stats, labels)
    if onto is None:
        onto = Ontology("label numbers only", ["labels"], {}, {})
    complete(onto, labels)
    onto.assign_ids()
    return onto


# --------------------------------------------------------------------------
# Per-region statistics
# --------------------------------------------------------------------------

# columns of the LightSuite tables that identify a row rather than measure it
ID_COLUMNS = {"sample", "atlas", "channel", "level", "region_id", "parcellation_index",
              "acronym", "name", "structure", "division", "hemisphere", "metric", "value"}


def attach_metrics(onto: Ontology, sample: Sample, label_count: dict) -> dict:
    """Adds numeric properties to onto.terms[key]['metrics'].

    Returns {property name: description}, in display order. Every numeric column of
    registration_analysis/region_metrics.csv is used (except ones that are constant,
    e.g. a channel's threshold), plus every metric of stats/region_stats_sample.csv
    (sample space) and stats/region_stats.csv (atlas space, prefixed atlas_).
    """
    for t in onto.terms.values():
        t["metrics"] = {}
    names = {}
    leaf_labels = {}
    for lab, chain in onto.chains.items():
        leaf_labels.setdefault(chain[0], []).append(lab)
    leaf_of = {labs[0]: k for k, labs in leaf_labels.items() if len(labs) == 1}

    rm_path = sample.root / "registration_analysis" / "region_metrics.csv"
    if rm_path.is_file():
        rm = pd.read_csv(rm_path)
        multi = rm["channel"].nunique() > 1
        cols = []
        for c in rm.columns:
            if c in ID_COLUMNS or not pd.api.types.is_numeric_dtype(rm[c]):
                continue
            if rm.groupby("channel")[c].nunique().max() <= 1:
                continue  # empty, or one value per channel (threshold, background)
            per_channel = multi and (rm.groupby(["level", "region_id"])[c].nunique() > 1).any()
            cols.append((c, per_channel))
        by_level_acr = {}
        for i, lv in enumerate(onto.levels):
            for chain in onto.chains.values():
                if chain[i] is not None:
                    by_level_acr[(lv, onto.terms[chain[i]]["acronym"])] = chain[i]
        for ch in sorted(rm["channel"].unique()):
            for c, per_channel in cols:
                if per_channel or ch == rm["channel"].min():
                    pname = f"{c}_ch{ch}" if per_channel else c
                    names[pname] = f"{c}{f', channel {ch}' if per_channel else ''} ({rm_path.name})"
        for _, r in rm.iterrows():
            if r["level"] == "label":
                key = leaf_of.get(int(r["region_id"]))
            else:
                key = by_level_acr.get((r["level"], r["acronym"]))
            if key is None:
                continue
            for c, per_channel in cols:
                if pd.notna(r[c]):
                    pname = f"{c}_ch{r['channel']}" if per_channel else c
                    onto.terms[key]["metrics"].setdefault(pname, float(r[c]))

    tables = [(sample.stats_sample, "", "sample space, stats/region_stats_sample.csv"),
              (sample.stats, "atlas_", "atlas space, stats/region_stats.csv")]
    for st, prefix, origin in tables:
        if st is None or not {"metric", "value", "hemisphere"} <= set(st):
            continue
        keycols = ["parcellation_index", "channel", "hemisphere", "metric"]
        if st.duplicated(keycols).any():
            continue  # per-segment tables (spinal cord)
        multi = st["channel"].nunique() > 1
        st = st[st["value"].notna()]
        for _, r in st.iterrows():
            key = leaf_of.get(int(r["parcellation_index"]))
            if key is None:
                continue
            m = str(r["metric"])
            per_channel = multi and not m.startswith("volume")
            ch = f"_ch{r['channel']}" if per_channel else ""
            pname = f"{prefix}{m}{ch}_{r['hemisphere']}"
            names.setdefault(pname, f"{m}, {r['hemisphere']} hemisphere"
                             f"{f', channel {r.channel}' if per_channel else ''} ({origin})")
            onto.terms[key]["metrics"][pname] = float(r["value"])
    used = {n for t in onto.terms.values() for n in t["metrics"]}
    return {n: d for n, d in names.items() if n in used}


# --------------------------------------------------------------------------
# Precomputed writers
# --------------------------------------------------------------------------

def write_precomputed(path: Path, arr_xyz: np.ndarray, seg: bool, voxel_nm, properties=False):
    path.mkdir(parents=True, exist_ok=True)
    dtype = "uint32" if seg else ("uint8" if arr_xyz.dtype == bool else str(arr_xyz.dtype))
    arr_xyz = arr_xyz.astype(dtype, copy=False)
    sx, sy, sz = arr_xyz.shape
    info = {
        "@type": "neuroglancer_multiscale_volume",
        "type": "segmentation" if seg else "image",
        "data_type": dtype,
        "num_channels": 1,
        "scales": [{
            "key": "s0", "size": [sx, sy, sz], "resolution": voxel_nm,
            "voxel_offset": [0, 0, 0], "chunk_sizes": [[CHUNK, CHUNK, CHUNK]],
            "encoding": "raw",
        }],
    }
    if properties:
        info["segment_properties"] = "segment_properties"
    (path / "info").write_text(json.dumps(info))
    sdir = path / "s0"
    sdir.mkdir()
    for x0 in range(0, sx, CHUNK):
        for y0 in range(0, sy, CHUNK):
            for z0 in range(0, sz, CHUNK):
                x1, y1, z1 = min(x0 + CHUNK, sx), min(y0 + CHUNK, sy), min(z0 + CHUNK, sz)
                # raw encoding: x fastest-varying (Fortran order)
                (sdir / f"{x0}-{x1}_{y0}-{y1}_{z0}-{z1}").write_bytes(
                    arr_xyz[x0:x1, y0:y1, z0:z1].tobytes(order="F"))


def tag_name(s: str) -> str:
    return re.sub(r"\W+", "_", s).strip("_") or "unassigned"


def write_segment_properties(path: Path, rows: list[dict], metric_names: dict):
    """rows: dicts with id, acronym, name, tags, metrics."""
    rows = sorted(rows, key=lambda r: r["id"])
    all_tags = sorted({t for r in rows for t in r["tags"]})
    tix = {t: k for k, t in enumerate(all_tags)}
    props = [
        {"id": "label", "type": "label", "values": [r["acronym"] for r in rows]},
        {"id": "description", "type": "description", "values": [r["name"] for r in rows]},
    ]
    if all_tags:
        props.append({"id": "tags", "type": "tags", "tags": all_tags,
                      "values": [sorted(tix[t] for t in r["tags"]) for r in rows]})
    for mname, desc in metric_names.items():
        # "NaN" (a string, JSON has no NaN) marks a missing value: Neuroglancer hides it
        props.append({"id": mname, "type": "number", "data_type": "float32",
                      "description": desc,
                      "values": [r["metrics"].get(mname, "NaN") for r in rows]})
    sp = path / "segment_properties"
    sp.mkdir(parents=True, exist_ok=True)
    (sp / "info").write_text(json.dumps({
        "@type": "neuroglancer_segment_properties",
        "inline": {"ids": [str(r["id"]) for r in rows], "properties": props},
    }))


def write_meshes(path: Path, annot: np.ndarray, term_labels: dict, voxel_nm) -> list[int]:
    """Legacy precomputed meshes: one fragment per term id, vertices in nm."""
    path.mkdir(parents=True)
    slices = ndimage.find_objects(annot.astype(np.int32, copy=False))
    n_lut = int(annot.max()) + 1
    vox = np.array(voxel_nm, dtype=np.float32)
    written = []
    for n, (tid, labels) in enumerate(term_labels.items(), 1):
        labels = [lab for lab in labels if lab < n_lut and slices[lab - 1] is not None]
        if not labels:
            continue
        boxes = [slices[lab - 1] for lab in labels]
        lo = [max(min(b[d].start for b in boxes) - 2, 0) for d in range(3)]
        hi = [min(max(b[d].stop for b in boxes) + 2, annot.shape[d]) for d in range(3)]
        crop = annot[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
        lut = np.zeros(n_lut, dtype=bool)
        lut[labels] = True
        mask = ndimage.gaussian_filter(np.pad(lut[crop], 1).astype(np.float32), 1.0)
        if mask.max() < 0.5:
            continue
        step = 1 if mask.size < 8e6 else 2
        verts, faces, _, _ = marching_cubes(mask, 0.5, step_size=step)
        # undo the 1-voxel pad; voxel centres sit at index + 0.5
        verts = (verts - 1 + np.array(lo) + 0.5) * vox
        # Neuroglancer asks for the manifest "<id>:0". Windows forbids ":" in file
        # names, so it is saved as "<id>_0" and CORSHandler maps the request.
        frag = f"{tid}_0_1"
        with open(path / frag, "wb") as f:
            f.write(struct.pack("<I", len(verts)))
            f.write(verts.astype("<f4").tobytes())
            f.write(faces.astype("<u4").tobytes())
        (path / f"{tid}_0").write_text(json.dumps({"fragments": [frag]}))
        written.append(int(tid))
        if n % 200 == 0:
            print(f"    {n}/{len(term_labels)}")
    (path / "info").write_text(json.dumps({
        "@type": "neuroglancer_legacy_mesh", "segment_properties": "segment_properties"}))
    return written


def region_centers(annot: np.ndarray, term_labels: dict) -> dict[str, list[float]]:
    """Per term id: a voxel inside the region near its centre (x, y, z, viewer voxel
    coordinates) followed by its bounding-box size in voxels. Used by the search panel
    to jump to a region."""
    slices = ndimage.find_objects(annot.astype(np.int32, copy=False))
    n_lut = int(annot.max()) + 1
    centers = {}
    for tid, labels in term_labels.items():
        labels = [lab for lab in labels if 0 < lab < n_lut and slices[lab - 1] is not None]
        if not labels:
            continue
        boxes = [slices[lab - 1] for lab in labels]
        lo = np.array([min(b[d].start for b in boxes) for d in range(3)])
        hi = np.array([max(b[d].stop for b in boxes) for d in range(3)])
        lut = np.zeros(n_lut, dtype=bool)
        lut[labels] = True
        mask = lut[annot[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]]
        c = np.array(ndimage.center_of_mass(mask))
        ci = tuple(np.clip(np.round(c).astype(int), 0, np.array(mask.shape) - 1))
        if not mask[ci]:  # non-convex region: take the region voxel nearest the centroid
            step = max(1, int(round((mask.size / 4e6) ** (1 / 3))))
            pts = np.argwhere(mask[::step, ::step, ::step]) * step
            if not len(pts):
                pts = np.argwhere(mask)
            c = pts[np.argmin(((pts - c) ** 2).sum(1))].astype(float)
        centers[str(tid)] = [round(float(v), 1) for v in (c + lo + 0.5)] + \
                            [int(v) for v in (hi - lo)]
    return centers


def extra_names(sample: Sample, onto: Ontology, key: str, values) -> dict[int, str]:
    names = {}
    if key == "segments" and sample.manifest.get("segments_csv"):
        p = Path(sample.manifest["segments_csv"])
        if p.is_file():
            segs = pd.read_csv(p)["Segment"].tolist()
            names = {i: segs[i - 1] for i in values if 0 < i <= len(segs)}
    if key == "division_labels":
        dirs = ([onto.dir] if onto.dir else []) + sample.atlas_dirs(None)
        legend = next((f for d in dirs for f in sorted(d.glob("division_id_legend*.csv"))), None)
        if legend is not None:
            df = pd.read_csv(legend)
            names = dict(zip(df["division_id"].astype(int), df["division_name"].astype(str)))
    return {i: names.get(i, f"{key} {i}") for i in values}


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------

def warn_label_overflow(onto: Ontology, annot: np.ndarray):
    """Atlas IDs that do not fit the annotation's integer type were corrupted before
    this tool saw them (e.g. Allen IDs like 312782574 passed through float32 and were
    then stored as uint16). The labels then match no region, and several regions share
    one label. This cannot be undone here, so say it clearly."""
    if not np.issubdtype(annot.dtype, np.integer):
        return
    top = np.iinfo(annot.dtype).max
    too_big = sorted(t["native_id"] for t in onto.terms.values()
                     if (t.get("native_id") or 0) > top and t["name"] and
                     not t["name"].endswith("(not in atlas)"))
    unknown = sorted(t["native_id"] for t in onto.terms.values()
                     if t["name"].endswith("(not in atlas)"))
    if too_big and unknown:
        print(f"  WARNING: the annotation volume is {annot.dtype} (max {top}), but the atlas has "
              f"{len(too_big)} region IDs above that (e.g. {too_big[0]}). Those regions were "
              f"corrupted when the annotation was saved, which is probably why {len(unknown)} labels "
              f"are 'not in atlas' ({', '.join(map(str, unknown[:6]))}"
              f"{', ...' if len(unknown) > 6 else ''}). Fix the annotation in the registration "
              "pipeline (keep it uint32, warp it with nearest-neighbour interpolation).")


def source_signature(sample: Sample, atlas_dir) -> str:
    h = hashlib.sha1(f"v{TOOL_VERSION}|{atlas_dir}".encode())
    h.update(Path(__file__).read_bytes())  # re-convert when this tool changes
    for f in sample.source_files():
        st = f.stat()
        h.update(f"{f}|{st.st_size}|{st.st_mtime_ns}".encode())
    return h.hexdigest()


def export(sample: Sample, out: Path, atlas_dir_override, meshes: bool) -> dict:
    print(f"Converting {sample.root} -> {out}")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    channels, images, annot, extra, voxel_nm = load_volumes(sample)
    print(f"  grid {annot.shape} (x,y,z), voxel {[v / 1000 for v in voxel_nm]} um, "
          f"channels {list(channels)}")

    meta = {"version": TOOL_VERSION, "sample": sample.name, "results": str(sample.root),
            "channels": list(channels), "ranges": {}, "levels": [], "extra": []}
    for name, arr in images.items():
        write_precomputed(out / name, arr, seg=False, voxel_nm=voxel_nm)
        nz = arr[::2, ::2, ::2]
        nz = nz[nz > 0] if (nz > 0).any() else nz.ravel()
        meta["ranges"][name] = [float(v) for v in np.percentile(nz, [1, 99.5])]
    meta["images"] = list(images)

    labels, counts = np.unique(annot, return_counts=True)
    label_count = dict(zip(labels.tolist(), counts.tolist()))
    onto = build_ontology(sample, labels, atlas_dir_override)
    print(f"  region names from: {onto.source}")
    warn_label_overflow(onto, annot)
    metric_names = attach_metrics(onto, sample, label_count)
    if metric_names:
        print(f"  region statistics: {', '.join(metric_names)}")
    else:
        print("  region statistics: none found")

    n_lut = max(int(annot.max()), max(onto.chains, default=0)) + 1
    term_labels, term_voxels = {}, {}
    level_terms = []
    for i, lv in enumerate(onto.levels):
        keys = {onto.chains[lab][i] for lab in onto.chains} - {None}
        level_terms.append(keys)
        for lab, chain in onto.chains.items():
            if chain[i] is None:
                continue
            tid = onto.ids[chain[i]]
            term_labels.setdefault(tid, set()).add(lab)
            if i == len(onto.levels) - 1:
                term_voxels[tid] = term_voxels.get(tid, 0) + label_count.get(lab, 0)
    # slice layers only for levels that actually split the tissue
    slice_levels = [i for i, keys in enumerate(level_terms) if len(keys) > 1]

    def tags_for(key, level_idx):
        chain = next(c for c in onto.chains.values() if c[level_idx] == key)
        tags = [tag_name(onto.terms[chain[j]]["acronym"]) for j in slice_levels
                if j >= max(level_idx + 1, 2) and chain[j] not in (None, key)]
        return sorted(set(tags + [onto.level_tags[level_idx]]))

    rows_by_id = {}
    for i, lv in enumerate(onto.levels):
        for key in level_terms[i]:
            tid = onto.ids[key]
            if tid not in rows_by_id:
                t = onto.terms[key]
                rows_by_id[tid] = {"id": tid, "acronym": t["acronym"], "name": t["name"],
                                   "tags": tags_for(key, i), "metrics": t["metrics"]}

    colors = {}
    for key, tid in onto.ids.items():
        if onto.terms[key].get("color"):
            colors[str(tid)] = onto.terms[key]["color"]

    for i in slice_levels:
        lv = onto.levels[i]
        lut = np.zeros(n_lut, dtype=np.uint32)
        for lab, chain in onto.chains.items():
            if chain[i] is not None:
                lut[lab] = onto.ids[chain[i]]
        vol = lut[annot]
        d = "level_" + onto.level_tags[i]
        write_precomputed(out / d, vol, seg=True, voxel_nm=voxel_nm, properties=True)
        ids = {onto.ids[k] for k in level_terms[i]}
        write_segment_properties(out / d, [rows_by_id[t] for t in ids], metric_names)
        meta["levels"].append({"name": lv, "dir": d})
        print(f"  level '{lv}': {len(ids)} regions")
    # region id at every level layer, per atlas label (for pinning a whole hierarchy)
    meta["hierarchy"] = sorted({tuple(onto.ids[c[i]] if c[i] is not None else 0
                                      for i in slice_levels) for c in onto.chains.values()})

    for key, vol in extra.items():
        if key == "division_labels" and "division" in onto.levels:
            continue  # same information as the division level
        vals = [int(v) for v in np.unique(vol) if v != 0]
        names = extra_names(sample, onto, key, vals)
        d = "extra_" + key
        write_precomputed(out / d, vol, seg=True, voxel_nm=voxel_nm, properties=True)
        write_segment_properties(out / d, [{"id": v, "acronym": names[v], "name": names[v],
                                            "tags": [], "metrics": {}} for v in vals], {})
        meta["extra"].append({"name": key, "dir": d})

    meta["mesh_default"] = []
    if meshes:
        print(f"  building region meshes ({len(term_labels)} regions) ...")
        written = write_meshes(out / "meshes", annot, term_labels, voxel_nm)
        write_segment_properties(out / "meshes", [rows_by_id[t] for t in written], metric_names)
        top = [t for t in written if t in {onto.ids[k] for k in level_terms[-1]}]
        if top:
            meta["mesh_default"] = [max(top, key=lambda t: term_voxels.get(t, 0))]
        print(f"  meshes: {len(written)}")
    meta["meshes"] = meshes
    meta["colors"] = colors
    print("  region centres (for search) ...")
    meta["centers"] = region_centers(annot, term_labels)
    return meta


# --------------------------------------------------------------------------
# Zoom: a LightSuite multires ROI (e.g. 3.2x) shown at native resolution
# --------------------------------------------------------------------------
#
# Nothing is registered or resampled here. The raw ROI TIFFs are converted as they
# are (s0 = native voxels) and placed with ONE affine built only from parameters the
# multires pipeline already saved:
#
#   ROI index --pair manifest--> ROI physical (um)
#             --roi_to_overview_tform (landmark fit)--> overview physical
#             --inverse of the elastix chain (TransformParameters.*.txt)--> overview physical
#             --pair manifest--> overview index (the 1x stack the atlas run used)
#             --registration-grid downsampling + permute_sample_to_atlas--> atlas viewer voxel
#
# The elastix chain maps fixed (overview) points to moving points, and the moving
# image is the ROI resampled through inv(roi_to_overview_tform), so
# overview_phys = inv(E) @ L @ roi_phys.

ZOOM_VERSION = 1
ZOOM_COLORS = {"488": "#00ff00", "561": "#ff00ff", "640": "#00ffff", "405": "#4080ff"}


def _resolve_near(p: str | None, *dirs: Path) -> Path | None:
    """A recorded path, or the same file name in one of *dirs* if the data moved."""
    if not p:
        return None
    p = Path(p)
    if p.exists():
        return p
    for d in dirs:
        if d is not None and (d / p.name).exists():
            return d / p.name
    return None


def _shift(v) -> np.ndarray:
    m = np.eye(4)
    m[:3, 3] = v
    return m


def _index_to_physical(spec: dict) -> np.ndarray:
    """lightsuite.multires.spec_geometry.index_xyz_to_physical as a 4x4 (x=column,
    y=row, z=page of the TIFF; index = voxel centre)."""
    m = np.eye(4)
    m[:3, :3] = np.asarray(spec["direction"], dtype=float).reshape(3, 3) @ np.diag(spec["spacing_um"])
    m[:3, 3] = spec["origin_um"]
    return m


def _read_elastix(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("(") and line.endswith(")"):
            key, *vals = line[1:-1].split()
            out[key] = [v.strip('"') for v in vals]
    return out


def _elastix_matrix(p: dict, path: Path) -> np.ndarray:
    """One elastix transform (fixed -> moving point) as a 4x4: A (x - c) + c + t."""
    kind = p["Transform"][0]
    par = [float(v) for v in p["TransformParameters"]]
    c = np.array([float(v) for v in p.get("CenterOfRotationPoint", [0, 0, 0])])
    if kind == "TranslationTransform":
        a, t = np.eye(3), np.array(par)
    elif kind == "EulerTransform":
        rx, ry, rz = par[:3]
        t = np.array(par[3:6])
        cx, sx, cy, sy, cz, sz = (np.cos(rx), np.sin(rx), np.cos(ry), np.sin(ry),
                                  np.cos(rz), np.sin(rz))
        Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
        Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
        zyx = p.get("ComputeZYX", ["false"])[0].lower() == "true"
        a = Rz @ Ry @ Rx if zyx else Rz @ Rx @ Ry  # itk::Euler3DTransform
    elif kind == "AffineTransform":
        a, t = np.array(par[:9]).reshape(3, 3), np.array(par[9:12])
    else:
        sys.exit(f"{path.name}: elastix '{kind}' is not affine, so the zoom cannot be "
                 "placed with one matrix (only Translation, Euler, Affine are supported).")
    m = np.eye(4)
    m[:3, :3] = a
    m[:3, 3] = t + c - a @ c
    return m


def _elastix_chain(last: Path) -> np.ndarray:
    """Composite fixed -> moving map of *last* and its InitialTransformParameterFileName chain."""
    p = _read_elastix(last)
    m = _elastix_matrix(p, last)
    init = p.get("InitialTransformParameterFileName", ["NoInitialTransform"])[0]
    if init != "NoInitialTransform":
        init_path = _resolve_near(init, last.parent)
        if init_path is None:
            sys.exit(f"{last.name}: initial transform {init} not found")
        how = p.get("HowToCombineTransforms", ["Compose"])[0]
        if how != "Compose":
            sys.exit(f"{last.name}: HowToCombineTransforms {how} is not supported")
        m = m @ _elastix_chain(init_path)  # T(x) = T_last(T_init(x))
    return m


def _registration_grid_matrix(sample: Sample, overview: dict) -> np.ndarray:
    """Overview-stack voxel (x=column, y=row, z=page) -> atlas viewer voxel (centres).

    Reproduces lightsuite.preprocess.slice_ops: XY = integer stride then a
    corner-aligned linear zoom to ceil(n * voxel / registres), Z = linear
    interpolation at native plane k * registres / voxel_z. Then the viewer's
    permute_brain_volume(perm) + to_xyz.
    """
    rp = sample.root / "regopts.json"
    if not rp.is_file():
        sys.exit(f"--zoom needs {rp} (native stack size and voxel size of the atlas run)")
    r = json.loads(rp.read_text(encoding="utf-8"))
    ny, nx, nz = int(r["ny"]), int(r["nx"]), int(r["nz"])
    vx, vy, vz = (float(v) for v in r["voxel_um"])
    res = float(r.get("registres_um", 20.0))

    o_nz, o_ny, o_nx = (int(v) for v in overview["shape_zyx"])
    if (o_nx, o_ny, o_nz) != (nx, ny, nz) or not np.allclose(overview["spacing_um"], [vx, vy, vz]):
        sys.exit(f"The zoom's overview stack {overview['volume_path']} "
                 f"({o_nx}x{o_ny}x{o_nz} at {overview['spacing_um']} um) is not the stack the "
                 f"atlas registration used ({nx}x{ny}x{nz} at {[vx, vy, vz]} um).")

    sxy, sz = vx / res, vz / res
    new_h, new_w, new_z = (max(1, int(np.ceil(ny * sxy))), max(1, int(np.ceil(nx * sxy))),
                           max(1, int(np.ceil(nz * sz))))

    def xy_factor(n, new):  # native index per grid index
        if np.isclose(sxy, 1.0) or new == n:
            return 1.0
        stride = max(1, min(ny // new_h, nx // new_w))
        n_s = len(range(0, n, stride))
        return float(stride) if n_s == new else stride * (n_s - 1) / (new - 1)

    fy, fx = xy_factor(ny, new_h), xy_factor(nx, new_w)
    fz = 1.0 / sz if new_z > 1 else 1.0
    # native (x, y, z) -> grid (Y, X, Z)
    g = np.zeros((4, 4))
    g[0, 1], g[1, 0], g[2, 2], g[3, 3] = 1 / fy, 1 / fx, 1 / fz, 1
    # permute_brain_volume on (Y, X, Z), then to_xyz swaps the first two axes
    n_grid = (new_h, new_w, new_z)
    p = np.zeros((4, 4))
    p[3, 3] = 1
    for d, v in enumerate(sample.perm):
        src = abs(v) - 1
        if v < 0:
            p[d, src], p[d, 3] = -1, n_grid[src] - 1
        else:
            p[d, src] = 1
    swap = np.eye(4)[[1, 0, 2, 3]]
    return swap @ p @ g


def _zoom_tag(label: str) -> str:
    tag = re.sub(r"(_\d+)+$", "", label or "")
    return re.sub(r"\W+", "_", tag) or "zoom"


class Zoom:
    """A multires results folder (multires_regopts.json) placed in the atlas viewer."""

    def __init__(self, root: Path, sample: Sample, atlas_info: dict):
        root = root.expanduser().resolve()
        found = [root] if (root / "multires_regopts.json").is_file() else \
            sorted(p.parent for p in root.glob("*/multires_regopts.json"))
        if len(found) != 1:
            sys.exit(f"--zoom {root}: expected one folder with multires_regopts.json, "
                     f"found {len(found)}")
        self.root = root = found[0]
        self.regopts = ro = json.loads((root / "multires_regopts.json").read_text(encoding="utf-8"))
        pair_path = _resolve_near(ro.get("pair_manifest_path"), root / "converted", root)
        if pair_path is None:
            sys.exit(f"{root}: pair manifest {ro.get('pair_manifest_path')} not found")
        self.pair = json.loads(pair_path.read_text(encoding="utf-8"))
        self.tag = _zoom_tag(ro.get("pair_label") or self.pair.get("pair_label", ""))

        # channels: name -> ROI spec (each with its own volume_path)
        chans = self.pair.get("channels") or {}
        self.channels = {str(k): dict(v["roi"]) for k, v in sorted(chans.items())}
        if not self.channels:
            self.channels = {str(self.pair.get("reference_channel") or "roi"): dict(self.pair["roi"])}
        for name, spec in self.channels.items():
            path = _resolve_near(spec["volume_path"])
            if path is None:
                sys.exit(f"{root}: ROI TIFF for channel {name} not found: {spec['volume_path']}")
            spec["_path"] = path

        # elastix chain (the last file holds the full InitialTransform chain)
        e_dir = _resolve_near(ro.get("elastix_output_dir"), root / "elastix_roi_to_overview")
        tps = [_resolve_near(t, e_dir) for t in ro.get("transform_paths") or []]
        if not tps and e_dir is not None:
            tps = sorted(e_dir.glob("TransformParameters.*.txt"))
        if not tps or any(t is None for t in tps):
            sys.exit(f"{root}: elastix TransformParameters files not found")
        self.transform_files = tps
        E = _elastix_chain(tps[-1])
        L = np.asarray(ro["roi_to_overview_tform"], float) if ro.get("roi_to_overview_tform") \
            else np.eye(4)  # geometry_mode 'metadata': stage positions only

        ref = str(self.pair.get("reference_channel") or next(iter(self.channels)))
        ref_pair = chans.get(ref, self.pair)
        overview, roi = ref_pair["overview"], ref_pair["roi"]
        for name, spec in self.channels.items():
            if not (np.allclose(spec["spacing_um"], roi["spacing_um"])
                    and np.allclose(spec["origin_um"], roi["origin_um"])
                    and list(spec["shape_zyx"]) == list(roi["shape_zyx"])):
                sys.exit(f"{root}: channel {name} ROI geometry differs from the reference "
                         f"channel {ref}; one transform cannot place both.")
        self.shape_zyx = [int(v) for v in roi["shape_zyx"]]

        # ROI voxel centre -> overview voxel centre -> atlas viewer voxel centre
        self.roi_to_overview_index = (np.linalg.inv(_index_to_physical(overview))
                                      @ np.linalg.inv(E) @ L @ _index_to_physical(roi))
        to_atlas = _registration_grid_matrix(sample, overview) @ self.roi_to_overview_index
        # Neuroglancer coordinates: voxel i spans [i, i+1), its centre is i + 0.5
        self.matrix = (_shift(0.5) @ to_atlas @ _shift(-0.5))[:3].tolist()
        self.atlas_scale_m = [v * 1e-9 for v in atlas_info["scales"][0]["resolution"]]
        nz_, ny_, nx_ = self.shape_zyx
        self.center = (np.array(self.matrix) @ np.array([nx_ / 2, ny_ / 2, nz_ / 2, 1.0])).tolist()
        corners = np.array([[x, y, z, 1.0] for x in (0, nx_) for y in (0, ny_) for z in (0, nz_)])
        pts = corners @ np.array(self.matrix).T
        self.extent = (pts.max(0) - pts.min(0)).tolist()  # ROI size in viewer voxels (x, y, z)

    def out_dir(self, out_arg: str | None) -> Path:
        if out_arg:
            return Path(out_arg).expanduser() / f"zoom_{self.tag}"
        # same rules as output_dir(): real write test, then the per-OS cache folder
        if _writable(self.root):
            return self.root / "neuroglancer"
        h = hashlib.sha1(str(self.root).encode()).hexdigest()[:8]
        return cache_root() / f"zoom_{self.tag}-{h}"

    def signature(self) -> str:
        h = hashlib.sha1(f"zoom{ZOOM_VERSION}|{CHUNK}".encode())
        for name, spec in self.channels.items():
            st = spec["_path"].stat()
            h.update(f"|{name}|{spec['_path']}|{st.st_size}|{st.st_mtime_ns}|"
                     f"{spec['spacing_um']}".encode())
        return h.hexdigest()


def _halve(v: np.ndarray) -> np.ndarray:
    """2x2x2 mean (Z, Y, X); odd sizes repeat the last plane (corner-aligned pyramid)."""
    Z, Y, X = v.shape
    out = np.empty(((Z + 1) // 2, (Y + 1) // 2, (X + 1) // 2), dtype=v.dtype)
    for z0 in range(0, Z, 64):
        b = v[z0:z0 + 64]
        b = np.pad(b, [(0, s % 2) for s in b.shape], mode="edge").astype(np.uint32)
        bz, by, bx = b.shape
        s = b.reshape(bz // 2, 2, by // 2, 2, bx // 2, 2).sum(axis=(1, 3, 5))
        out[z0 // 2:z0 // 2 + bz // 2] = ((s + 4) // 8).astype(v.dtype)
    return out


def _load_stack_zyx(path: Path) -> np.ndarray:
    try:
        return np.array(tifffile.memmap(path, mode="r"))
    except (ValueError, OSError):
        with tifffile.TiffFile(path) as tf:
            return np.stack([p.asarray() for p in tf.pages])


def write_zoom_channel(path: Path, tif: Path, voxel_um_xyz) -> list[float]:
    """Raw multi-page TIFF -> precomputed multiscale image (s0 = native voxels)."""
    from concurrent.futures import ThreadPoolExecutor

    print(f"  reading {tif.name} ...")
    vol = _load_stack_zyx(tif)
    if vol.ndim != 3:
        sys.exit(f"{tif}: expected a 3D stack, got shape {vol.shape}")
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    scales, level, k, sample_vals = [], vol, 0, None
    with ThreadPoolExecutor(max_workers=min(16, os.cpu_count() or 4)) as pool:
        while True:
            Z, Y, X = level.shape
            key = f"s{k}"
            (path / key).mkdir()

            def write(z0, level=level, key=key, Z=Z, Y=Y, X=X):
                for y0 in range(0, Y, CHUNK):
                    for x0 in range(0, X, CHUNK):
                        z1, y1, x1 = min(z0 + CHUNK, Z), min(y0 + CHUNK, Y), min(x0 + CHUNK, X)
                        # (Z, Y, X) C order == x fastest-varying, as raw precomputed wants
                        (path / key / f"{x0}-{x1}_{y0}-{y1}_{z0}-{z1}").write_bytes(
                            np.ascontiguousarray(level[z0:z1, y0:y1, x0:x1]).tobytes())

            list(pool.map(write, range(0, Z, CHUNK)))
            f = 2 ** k
            scales.append({"key": key, "size": [X, Y, Z],
                           "resolution": [float(v) * 1000 * f for v in voxel_um_xyz],
                           "voxel_offset": [0, 0, 0], "chunk_sizes": [[CHUNK, CHUNK, CHUNK]],
                           "encoding": "raw"})
            print(f"    {key}: {X}x{Y}x{Z}")
            if sample_vals is None and (k == 2 or max(X, Y, Z) <= 2 * CHUNK):
                sample_vals = level[level > 0]
            if max(X, Y, Z) <= 2 * CHUNK:
                break
            level, k = _halve(level), k + 1
    info = {"@type": "neuroglancer_multiscale_volume", "type": "image",
            "data_type": str(vol.dtype), "num_channels": 1, "scales": scales}
    (path / "info").write_text(json.dumps(info))  # written last: marks a complete conversion
    if sample_vals is None or not sample_vals.size:
        sample_vals = level.ravel()
    return [float(v) for v in np.percentile(sample_vals, [1, 99.9])]


def prepare_zoom(zoom: Zoom, out_arg: str | None, reexport: bool) -> dict:
    out = zoom.out_dir(out_arg)
    meta_path = out / "zoom_meta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
    sig = zoom.signature()
    if reexport or meta.get("signature") != sig or \
            not all((out / ch / "info").is_file() for ch in zoom.channels):
        print(f"Converting zoom {zoom.root} -> {out}")
        out.mkdir(parents=True, exist_ok=True)
        meta = {"version": ZOOM_VERSION, "ranges": {}}
        for ch, spec in zoom.channels.items():
            meta["ranges"][ch] = write_zoom_channel(out / ch, spec["_path"], spec["spacing_um"])
        meta["signature"] = sig
    else:
        print(f"Using cached zoom conversion in {out}")
    # the placement is cheap and always recomputed from the saved parameters
    meta.update({"root": str(zoom.root), "tag": zoom.tag, "channels": list(zoom.channels),
                 "matrix_roi_voxel_to_atlas_viewer_voxel": zoom.matrix,
                 "atlas_viewer_voxel_m": zoom.atlas_scale_m,
                 "transform_files": [str(t) for t in zoom.transform_files]})
    meta_path.write_text(json.dumps(meta, indent=1))
    meta["dir"] = out
    return meta


# --------------------------------------------------------------------------
# Viewer
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Channel look-up tables (LUTs)
# --------------------------------------------------------------------------
#
# Each LUT is turned into a Neuroglancer shader. The "normalized" invlerp control
# (contrast range) is kept when the LUT changes. Colormaps are 33 samples of the
# matplotlib maps, linearly interpolated in the shader (error < 1%).

_TABLES = {  # 33 samples of the matplotlib colormaps, as 0-255 RGB
    "inferno": ("000004 040312 0b0724 150b37 210c4a 2f0a5b 3d0965 4a0c6b 57106e 64156e 71196e "
             "7d1e6d 8a226a 972766 a32c61 b0315b bc3754 c73e4c d24644 db503b e45a31 eb6628 "
             "f1731d f68013 f98e09 fb9d07 fcac11 fbbc21 f9cb35 f5db4c f2ea69 f3f68a fcffa4").split(),
    "magma": ("000004 030312 0a0822 130d34 1d1147 29115a 36106b 440f76 51127c 5d177f 6a1c81 "
             "762181 832681 902a81 9c2e7f aa337d b73779 c43c75 d0416f dc4869 e75263 ef5d5e "
             "f56b5c f9795d fc8961 fd9869 fea772 feb67c fec488 fed395 fde2a3 fcf0b2 fcfdbf").split(),
    "plasma": ("0d0887 220690 310597 3f049c 4c02a1 5901a5 6600a7 7201a8 7e03a8 8a09a5 9511a1 "
             "a01a9c aa2395 b32c8e bc3587 c43e7f cc4778 d35171 da5a6a e06363 e66c5c eb7655 "
             "f0804e f58b47 f89540 fba139 fdac33 feb82c fdc527 fcd225 f8df25 f4ed27 f0f921").split(),
    "viridis": ("440154 470d60 48186a 482374 472d7b 453781 424086 3e4989 3b528b 375b8d 33638d "
             "2f6b8e 2c728e 297a8e 26828e 23898e 21918c 1f988b 1fa088 22a785 28ae80 32b67a "
             "3fbc73 4ec36b 5ec962 70cf57 84d44b 98d83e addc30 c2df23 d8e219 ece51b fde725").split(),
    "turbo": ("30123b 392a73 4040a2 4456c7 466be3 4680f6 4294ff 37a8fa 28bceb 1ccdd8 18ddc2 "
             "1fe9af 32f298 4ef97d 6dfe62 8bff4b a4fc3c b9f635 cdec34 dfdf37 eecf3a f8be39 "
             "fdac34 fe962b fb7e21 f46617 eb500e df3f08 d02f05 be2102 a91601 920b01 7a0403").split(),
}


# name -> (label, kind, parameter). kind: "color" (black -> colour), "table", "fn"
LUTS = {
    "grays": ("Grays", "color", "#ffffff"),
    "green": ("Green", "color", "#00ff00"),
    "red": ("Red", "color", "#ff0000"),
    "blue": ("Blue", "color", "#0060ff"),
    "magenta": ("Magenta", "color", "#ff00ff"),
    "cyan": ("Cyan", "color", "#00ffff"),
    "yellow": ("Yellow", "color", "#ffff00"),
    "fire": ("Fire (hot)", "fn", "fire"),
    "inferno": ("Inferno", "table", "inferno"),
    "magma": ("Magma", "table", "magma"),
    "plasma": ("Plasma", "table", "plasma"),
    "viridis": ("Viridis", "table", "viridis"),
    "turbo": ("Turbo", "table", "turbo"),
    "jet": ("Jet", "fn", "jet"),
    "hilo": ("HiLo (saturation check)", "fn", "hilo"),
}
_LUT_FOR_COLOR = {"#ffffff": "grays", "#00ff00": "green", "#ff00ff": "magenta",
                  "#00ffff": "cyan", "#ffff00": "yellow", "#4080ff": "blue", "#ff0000": "red"}


def _hex_rgb(h: str) -> np.ndarray:
    return np.array([int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def lut_rgb(lut: str, t: np.ndarray) -> np.ndarray:
    """The LUT evaluated in numpy (same formulas as the shader), for the panel swatches."""
    _, kind, par = LUTS[lut]
    t = np.clip(np.asarray(t, dtype=float), 0, 1)[:, None]
    if kind == "color":
        rgb = t * _hex_rgb(par)
    elif kind == "table":
        tab = np.array([_hex_rgb("#" + h) for h in _TABLES[par]])
        x = np.linspace(0, 1, len(tab))
        rgb = np.stack([np.interp(t[:, 0], x, tab[:, k]) for k in range(3)], axis=1)
    elif par == "fire":
        rgb = np.hstack([3 * t, 3 * t - 1, 3 * t - 2])
    elif par == "jet":
        rgb = np.hstack([1.5 - np.abs(4 * t - 3), 1.5 - np.abs(4 * t - 2), 1.5 - np.abs(4 * t - 1)])
    else:  # hilo: gray, with blue at the bottom and red at the top of the range
        rgb = np.repeat(t, 3, axis=1)
        rgb[t[:, 0] <= 0.0] = (0, 0, 1)
        rgb[t[:, 0] >= 1.0] = (1, 0, 0)
    return np.clip(rgb, 0, 1)


def lut_css(lut: str, invert: bool = False) -> str:
    t = np.linspace(0, 1, 9)
    if lut == "hilo":  # show the clipping colours as narrow end bands
        t = np.array([0, 0.001, 0.25, 0.5, 0.75, 0.999, 1.0])
    rgb = lut_rgb(lut, 1 - t if invert else t)
    stops = ", ".join(f"rgb({int(r * 255)},{int(g * 255)},{int(b * 255)}) {p * 100:.1f}%"
                      for (r, g, b), p in zip(rgb, t))
    return f"linear-gradient(90deg, {stops})"


def lut_shader(lut: str, invert: bool = False) -> str:
    _, kind, par = LUTS[lut]
    head = ""
    if kind == "color":
        head = f'#uicontrol vec3 color color(default="{par}")\n'
        body = "vec3 rgb = color * t;"
    elif kind == "table":
        tab = _TABLES[par]
        vals = ", ".join("vec3(%.4f, %.4f, %.4f)" % tuple(_hex_rgb("#" + h)) for h in tab)
        n = len(tab) - 1
        head = f"const vec3 LUT[{n + 1}] = vec3[{n + 1}]({vals});\n"
        body = (f"float x = t * {n}.0; int i = int(min(floor(x), {n - 1}.0));"
                " vec3 rgb = mix(LUT[i], LUT[i + 1], x - float(i));")
    elif par == "fire":
        body = "vec3 rgb = clamp(vec3(3.0 * t, 3.0 * t - 1.0, 3.0 * t - 2.0), 0.0, 1.0);"
    elif par == "jet":
        body = ("vec3 rgb = clamp(vec3(1.5 - abs(4.0 * t - 3.0), 1.5 - abs(4.0 * t - 2.0),"
                " 1.5 - abs(4.0 * t - 1.0)), 0.0, 1.0);")
    else:
        body = ("vec3 rgb = vec3(t); if (t <= 0.0) rgb = vec3(0.0, 0.0, 1.0);"
                " if (t >= 1.0) rgb = vec3(1.0, 0.0, 0.0);")
    t = "1.0 - clamp(normalized(), 0.0, 1.0)" if invert else "clamp(normalized(), 0.0, 1.0)"
    return (f"// lightsuite-lut: {lut}{' inverted' if invert else ''}\n{head}"
            "#uicontrol invlerp normalized\n"
            f"void main() {{\n  float t = {t};\n  {body}\n  emitRGB(rgb);\n}}\n")


# --------------------------------------------------------------------------
# Region search
# --------------------------------------------------------------------------

def _norm(s: str) -> str:
    """Lower case, no accents/punctuation, British spellings folded (coeruleus = ceruleus)."""
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = s.replace("ae", "e").replace("oe", "e")
    return " ".join(s.split())


class RegionIndex:
    """Every region of every hierarchy level / label layer, searchable by name or acronym."""

    def __init__(self, out: Path, meta: dict):
        self.entries = []
        seen = set()
        layers = [(lv, True) for lv in meta["levels"]] + [(ex, False) for ex in meta["extra"]]
        for rank, (lv, is_level) in enumerate(layers):
            info = out / lv["dir"] / "segment_properties" / "info"
            if not info.is_file():
                continue
            inline = json.loads(info.read_text(encoding="utf-8"))["inline"]
            props = {p["id"]: p["values"] for p in inline["properties"]
                     if p["id"] in ("label", "description")}
            labels = props.get("label", [""] * len(inline["ids"]))
            names = props.get("description", labels)
            for sid, acr, name in zip(inline["ids"], labels, names):
                key = (acr.lower(), name.lower())
                if key in seen:  # same region at a coarser level: keep the finest one
                    continue
                seen.add(key)
                self.entries.append({
                    "layer": lv["name"], "id": int(sid), "acronym": acr, "name": name,
                    "rank": rank, "is_level": is_level,
                    "_n": _norm(name), "_a": _norm(acr), "_toks": _norm(name).split()})

    @staticmethod
    def _score(q: str, qtoks: list[str], e: dict) -> float:
        n, a = e["_n"], e["_a"]
        if q == a or q == n:
            return 100
        if n.startswith(q):
            return 92
        if a.startswith(q) and len(q) >= 2:
            return 88
        if f" {q}" in f" {n}":
            return 84
        if q in n:
            return 78
        toks = e["_toks"]
        if all(any(t.startswith(qt) for t in toks) for qt in qtoks):
            return 74
        # fuzzy: every typed word close to some word of the name (typos, spellings)
        sims = [max((difflib.SequenceMatcher(None, qt, t).ratio() for t in toks), default=0)
                for qt in qtoks]
        fuzzy = sum(sims) / len(sims) if sims else 0
        whole = difflib.SequenceMatcher(None, q, n).ratio()
        best = max(fuzzy if min(sims, default=0) >= 0.6 else 0, whole)
        return 70 * best if best >= 0.72 else 0

    def search(self, query: str, limit: int = 25) -> list[dict]:
        q = _norm(query)
        if not q:
            return []
        if q.isdigit():
            hits = [(100, e) for e in self.entries if e["id"] == int(q)]
        else:
            qtoks = q.split()
            hits = [(s, e) for e in self.entries if (s := self._score(q, qtoks, e)) > 0]
        hits.sort(key=lambda se: (-se[0], se[1]["rank"], len(se[1]["name"])))
        return [{k: v for k, v in e.items() if not k.startswith("_")} | {"score": round(s)}
                for s, e in hits[:limit]]


# --------------------------------------------------------------------------
# Control panel (search + channel LUTs) next to the Neuroglancer viewer
# --------------------------------------------------------------------------

class ViewerApp:
    """Server-side actions of the panel; all of them edit the Neuroglancer viewer state."""

    def __init__(self, viewer, meta: dict, channel_layers: list[str], level_names: list[str],
                 seg_names: list[str]):
        self.viewer, self.meta = viewer, meta
        self.channels = channel_layers
        self.level_names, self.seg_names = level_names, seg_names
        self.token = secrets.token_urlsafe(16)
        self.luts = {}  # layer -> {"lut", "invert"}
        self.highlight = []  # [(layer, id, label)]
        self._saved_vis, self._saved_alpha, self._saved_view = None, {}, (None, None)
        self.index: RegionIndex | None = None
        self.viewer_url = ""
        self.lock = threading.Lock()

    # --- hierarchy -------------------------------------------------------
    def pins_for(self, layer: str, seg: int) -> dict:
        """The region and the regions containing it, at every level layer."""
        pinned = {layer: seg}
        if layer in self.level_names:
            levels = [self.level_names.index(layer)]
        elif layer == "region meshes":
            levels = range(len(self.level_names))
        else:
            return pinned
        for k in levels:
            row = next((r for r in self.meta["hierarchy"] if r[k] == seg), None)
            if row is not None:
                for j in range(k, len(self.level_names)):
                    if row[j]:
                        pinned.setdefault(self.level_names[j], row[j])
                break
        return pinned

    def _selection(self, pinned: dict, position) -> dict:
        order = {n: i for i, n in enumerate(self.seg_names)}
        return {"position": [float(x) for x in position],
                "layers": {n: {"value": str(pinned[n])}
                           for n in sorted(pinned, key=lambda n: order.get(n, 99))}}

    def _status(self, msg: str | None):
        with self.viewer.config_state.txn() as c:
            if msg:
                c.status_messages["lightsuite"] = msg
            else:
                c.status_messages.pop("lightsuite", None)

    # --- search ------------------------------------------------------------
    def show_region(self, layer: str, seg: int, add: bool = False) -> dict:
        center = self.meta.get("centers", {}).get(str(seg))
        with self.lock, self.viewer.txn() as s:
            if layer not in s.layers:
                return {"error": f"no layer {layer}"}
            seg_layers = [n for n in self.seg_names if n in s.layers and n != "region meshes"]
            if self._saved_vis is None:  # view before the first search, restored by clear()
                self._saved_vis = {n: s.layers[n].visible for n in seg_layers}
                self._saved_alpha = {n: s.layers[n].selected_alpha for n in seg_layers}
                self._saved_view = (s.position, s.cross_section_scale)
            if not add:
                self.highlight = []
                for n in seg_layers:
                    s.layers[n].segments = set()
                    s.layers[n].visible = False
            self.highlight = [h for h in self.highlight if h[:2] != (layer, seg)]
            label = next((e for e in self.index.entries if e["layer"] == layer and e["id"] == seg),
                         {"name": str(seg), "acronym": ""})
            self.highlight.append((layer, seg, f"{label['name']} ({label['acronym']})"))
            for n, sid, _ in self.highlight:
                s.layers[n].visible = True
                s.layers[n].selected_alpha = 0.7
                s.layers[n].segments.add(sid)
            if "region meshes" in s.layers:
                s.layers["region meshes"].segments = {
                    sid for n, sid, _ in self.highlight if n in self.level_names}
            if center:
                s.position = center[:3]
                extent = max(center[3:6])
                s.cross_section_scale = float(np.clip(extent * 1.8 / 500, 0.12, 1.5))
            pinned = {}
            for n, sid, _ in self.highlight:
                for k, v in self.pins_for(n, sid).items():
                    pinned.setdefault(k, v)
            s._json_data["selection"] = self._selection(pinned, center[:3] if center else
                                                        s.position)
            s.selected_layer.layer = layer
        names = "; ".join(h[2] for h in self.highlight)
        self._status(f"Highlighted: {names}")
        return {"ok": True, "highlight": [h[2] for h in self.highlight], "found": bool(center)}

    def clear(self) -> dict:
        with self.lock, self.viewer.txn() as s:
            for n, vis in (self._saved_vis or {}).items():
                if n in s.layers:
                    s.layers[n].segments = set()
                    s.layers[n].visible = vis
                    s.layers[n].selected_alpha = self._saved_alpha[n]
            if "region meshes" in s.layers:
                s.layers["region meshes"].segments = set(self.meta.get("mesh_default", []))
            if self._saved_vis is not None and self._saved_view[0] is not None:
                s.position, s.cross_section_scale = self._saved_view
            s._json_data["selection"] = {}
            self._saved_vis, self.highlight = None, []
        self._status(None)
        return {"ok": True}

    # --- channels ----------------------------------------------------------
    def set_lut(self, layer: str, lut: str, invert: bool) -> dict:
        if layer not in self.channels or lut not in LUTS:
            return {"error": "unknown layer or LUT"}
        with self.lock, self.viewer.txn() as s:
            L = s.layers[layer]
            ctl = L.to_json().get("shaderControls", {})
            new = {"normalized": ctl["normalized"]} if "normalized" in ctl else {}
            if LUTS[lut][1] == "color":
                new["color"] = LUTS[lut][2]
            L.shader = lut_shader(lut, invert)
            L.shader_controls = new
        self.luts[layer] = {"lut": lut, "invert": invert}
        return {"ok": True}

    def set_visible(self, layer: str, visible: bool) -> dict:
        if layer not in self.channels:
            return {"error": "unknown layer"}
        with self.lock, self.viewer.txn() as s:
            s.layers[layer].visible = visible
        return {"ok": True}

    def state(self) -> dict:
        st = self.viewer.state
        chans = []
        for n in self.channels:
            cur = self.luts.get(n, {"lut": "grays", "invert": False})
            chans.append({"layer": n, "visible": bool(st.layers[n].visible) if n in st.layers
                          else False, **cur})
        return {"channels": chans, "highlight": [h[2] for h in self.highlight]}

    def config(self, viewer_url: str) -> dict:
        return {"token": self.token, "viewer_url": viewer_url, "sample": self.meta["sample"],
                "luts": [{"id": k, "label": v[0], "css": lut_css(k),
                          "css_inv": lut_css(k, True)} for k, v in LUTS.items()],
                "regions": len(self.index.entries)}

    # --- HTTP ------------------------------------------------------------------
    def handle(self, method: str, path: str, query: dict, body: dict):
        """Returns (status, content_type, bytes) for /lightsuite/... requests."""
        if path in ("/lightsuite", "/lightsuite/"):
            html = PANEL_HTML.replace("/*CONFIG*/null", json.dumps(self.config(self.viewer_url)))
            return 200, "text/html; charset=utf-8", html.encode("utf-8")
        tok = query.get("token", [""])[0] or body.get("token", "")
        if not secrets.compare_digest(str(tok), self.token):
            return 403, "application/json", b'{"error": "bad token"}'
        api = path.removeprefix("/lightsuite/api/")
        if api == "search":
            res = {"results": self.index.search(query.get("q", [""])[0])}
        elif api == "state":
            res = self.state()
        elif method == "POST" and api == "show":
            res = self.show_region(str(body["layer"]), int(body["id"]), bool(body.get("add")))
        elif method == "POST" and api == "clear":
            res = self.clear()
        elif method == "POST" and api == "lut":
            res = self.set_lut(str(body["layer"]), str(body["lut"]), bool(body.get("invert")))
        elif method == "POST" and api == "visible":
            res = self.set_visible(str(body["layer"]), bool(body["visible"]))
        else:
            return 404, "application/json", b'{"error": "not found"}'
        return 200, "application/json", json.dumps(res).encode("utf-8")


PANEL_HTML = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>LightSuite Viewer</title>
<style>
:root { --bg:#111315; --panel:#1a1d20; --line:#2c3034; --fg:#e6e7e8; --muted:#9aa0a6;
        --accent:#6cb8e0; --hit:#23303a; }
* { box-sizing:border-box; }
html, body { margin:0; height:100%; background:var(--bg); color:var(--fg);
             font:14px/1.4 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
.app { display:flex; height:100vh; }
aside { width:310px; flex:none; background:var(--panel); border-right:1px solid var(--line);
        display:flex; flex-direction:column; overflow:hidden; transition:width .15s; }
.app.collapsed aside { width:0; border:0; }
.inner { width:310px; padding:12px 14px; overflow-y:auto; flex:1; }
iframe { flex:1; border:0; min-width:0; background:#000; }
h1 { font-size:15px; margin:0 0 2px; }
.sample { color:var(--muted); font-size:12px; margin-bottom:12px; word-break:break-all; }
h2 { font-size:12px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted);
     margin:16px 0 6px; font-weight:600; }
input[type=search] { width:100%; padding:8px 10px; border-radius:6px; border:1px solid var(--line);
     background:#0f1113; color:var(--fg); font-size:14px; }
input[type=search]:focus { outline:2px solid var(--accent); outline-offset:-1px; }
ul.results { list-style:none; margin:6px 0 0; padding:0; max-height:42vh; overflow-y:auto; }
ul.results li { padding:6px 8px; border-radius:5px; cursor:pointer; }
ul.results li.active, ul.results li:hover { background:var(--hit); }
.rname { display:block; }
.rmeta { color:var(--muted); font-size:12px; }
.hint { color:var(--muted); font-size:12px; margin:6px 0 0; }
.current { margin-top:8px; padding:8px; border-radius:6px; background:#162229;
           border:1px solid #24404f; font-size:13px; display:none; }
.current.on { display:block; }
button { background:#2a2f34; color:var(--fg); border:1px solid var(--line); border-radius:5px;
         padding:4px 10px; cursor:pointer; font-size:13px; }
button:hover { background:#343a40; }
.chan { padding:8px 0; border-bottom:1px solid var(--line); }
.chan-top { display:flex; align-items:center; gap:8px; }
.chan-top label.name { flex:1; font-weight:600; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
select { flex:1; background:#0f1113; color:var(--fg); border:1px solid var(--line); border-radius:5px;
         padding:4px 6px; font-size:13px; min-width:0; }
.chan-row { display:flex; gap:8px; align-items:center; margin-top:6px; }
.swatch { height:10px; border-radius:3px; margin-top:6px; border:1px solid var(--line); }
.inv { color:var(--muted); font-size:12px; white-space:nowrap; }
.toggle { position:fixed; left:0; top:50%; z-index:5; transform:translateY(-50%);
          border-radius:0 6px 6px 0; padding:10px 4px; display:none; }
.app.collapsed .toggle { display:block; }
.err { color:#f28b82; }
@media (max-width: 700px) {
  .app { flex-direction:column; }
  aside { width:100%; max-height:45vh; border-right:0; border-bottom:1px solid var(--line); }
  .inner { width:100%; }
  .app.collapsed aside { width:100%; max-height:0; }
  .toggle { top:0; left:50%; transform:translateX(-50%); border-radius:0 0 6px 6px; padding:2px 12px; }
}
</style></head>
<body>
<div class="app" id="app">
  <aside><div class="inner">
    <div style="display:flex;justify-content:space-between;align-items:start;gap:8px">
      <div><h1>LightSuite viewer</h1><div class="sample" id="sample"></div></div>
      <button id="hide" title="Hide panel">&#x2039;</button>
    </div>

    <h2>Find a region</h2>
    <input type="search" id="q" placeholder="e.g. locus coeruleus, CA1, Isocortex" autocomplete="off" autofocus>
    <ul class="results" id="results"></ul>
    <p class="hint" id="hint">Type a name, acronym or ID. Enter shows the first match,
      &uarr;&darr; to choose. Shift+click adds a region.</p>
    <div class="current" id="current">
      <div id="curtext"></div>
      <div style="margin-top:6px"><button id="clear">Show all regions</button></div>
    </div>

    <h2>Channels</h2>
    <div id="chans"></div>
    <p class="hint">Contrast: use the range slider in the layer's Rendering tab in Neuroglancer.</p>
  </div></aside>
  <button class="toggle" id="show" title="Show panel">&#x203A;</button>
  <iframe id="ng" title="Neuroglancer" allow="fullscreen"></iframe>
</div>
<script>
const CFG = /*CONFIG*/null;
const $ = id => document.getElementById(id);
$("sample").textContent = CFG.sample + " · " + CFG.regions + " regions";
$("ng").src = CFG.viewer_url;
const api = (p, body) => fetch("/lightsuite/api/" + p + (body ? "" : (p.includes("?") ? "&" : "?") + "token=" + CFG.token),
  body ? {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({...body, token:CFG.token})} : {})
  .then(r => r.json());

// ---- search
let results = [], active = -1, timer = null, seq = 0;
function render() {
  const ul = $("results"); ul.innerHTML = "";
  results.forEach((r, i) => {
    const li = document.createElement("li");
    if (i === active) li.className = "active";
    const a = document.createElement("span"); a.className = "rname"; a.textContent = r.name;
    const m = document.createElement("span"); m.className = "rmeta";
    m.textContent = (r.acronym && r.acronym !== r.name ? r.acronym + " · " : "") + r.layer +
                    (r.score < 70 ? " · similar" : "");
    li.append(a, m);
    li.onclick = ev => choose(i, ev.shiftKey);
    ul.append(li);
  });
}
function doSearch() {
  const q = $("q").value.trim(), my = ++seq;
  if (!q) { results = []; active = -1; render(); $("hint").classList.remove("err"); return; }
  api("search?q=" + encodeURIComponent(q)).then(d => {
    if (my !== seq) return;
    results = d.results || []; active = results.length ? 0 : -1; render();
    $("hint").textContent = results.length ? results.length + (results.length === 25 ? "+" : "") + " matches"
                                            : "No region matches “" + q + "”.";
    $("hint").classList.toggle("err", !results.length);
  });
}
function choose(i, add) {
  const r = results[i]; if (!r) return;
  active = i; render();
  api("show", {layer:r.layer, id:r.id, add:!!add}).then(showCurrent);
}
function showCurrent(d) {
  const on = d && d.highlight && d.highlight.length;
  $("current").classList.toggle("on", !!on);
  if (on) $("curtext").textContent = "Highlighted: " + d.highlight.join("; ");
}
$("q").addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(doSearch, 120); });
$("q").addEventListener("keydown", ev => {
  if (ev.key === "ArrowDown") { active = Math.min(active + 1, results.length - 1); render(); ev.preventDefault(); }
  else if (ev.key === "ArrowUp") { active = Math.max(active - 1, 0); render(); ev.preventDefault(); }
  else if (ev.key === "Enter") { clearTimeout(timer); if (results.length) choose(Math.max(active, 0), ev.shiftKey); else doSearch(); }
  else if (ev.key === "Escape") { $("q").value = ""; doSearch(); }
});
$("clear").onclick = () => api("clear", {}).then(() => showCurrent(null));

// ---- channels
const LUTS = Object.fromEntries(CFG.luts.map(l => [l.id, l]));
function buildChannels(st) {
  const box = $("chans"); box.innerHTML = "";
  st.channels.forEach(c => {
    const div = document.createElement("div"); div.className = "chan";
    const id = "v_" + c.layer.replace(/\W/g, "_");
    div.innerHTML = `<div class="chan-top"><input type="checkbox" id="${id}">
        <label class="name" for="${id}"></label></div>
      <div class="chan-row"><select aria-label="LUT"></select>
        <label class="inv"><input type="checkbox"> invert</label></div>
      <div class="swatch"></div>`;
    div.querySelector("label.name").textContent = c.layer;
    const vis = div.querySelector("#" + id), sel = div.querySelector("select"),
          inv = div.querySelector(".inv input"), sw = div.querySelector(".swatch");
    CFG.luts.forEach(l => sel.add(new Option(l.label, l.id)));
    vis.checked = c.visible; sel.value = c.lut; inv.checked = c.invert;
    const paint = () => sw.style.background = LUTS[sel.value][inv.checked ? "css_inv" : "css"];
    paint();
    const send = () => { paint(); api("lut", {layer:c.layer, lut:sel.value, invert:inv.checked}); };
    sel.onchange = send; inv.onchange = send;
    vis.onchange = () => api("visible", {layer:c.layer, visible:vis.checked});
    div.dataset.layer = c.layer;
    box.append(div);
  });
}
function syncVisibility() {
  api("state").then(st => st.channels.forEach(c => {
    const d = document.querySelector(`.chan[data-layer="${CSS.escape(c.layer)}"] input[type=checkbox]`);
    if (d) d.checked = c.visible;
  })).catch(() => {});
}
api("state").then(st => { buildChannels(st); showCurrent(st); });
setInterval(syncVisibility, 3000);

$("hide").onclick = () => $("app").classList.add("collapsed");
$("show").onclick = () => $("app").classList.remove("collapsed");
</script>
</body></html>
"""


class CORSHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, mounts=None, app=None, **kwargs):
        self.mounts = mounts or {}  # "/zoom/<tag>" -> folder served under that prefix
        self.app = app  # ViewerApp: the /lightsuite/ panel and its API
        super().__init__(*args, **kwargs)

    def _app_request(self, method: str) -> bool:
        url = urllib.parse.urlsplit(self.path)
        if self.app is None or not url.path.startswith("/lightsuite"):
            return False
        body = {}
        if method == "POST":
            n = int(self.headers.get("Content-Length") or 0)
            try:
                body = json.loads(self.rfile.read(n) or b"{}")
            except ValueError:
                body = {}
        try:
            code, ctype, data = self.app.handle(method, url.path, urllib.parse.parse_qs(url.query),
                                                body)
        except (KeyError, ValueError, TypeError) as e:
            code, ctype, data = 400, "application/json", json.dumps({"error": str(e)}).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)
        return True

    def do_GET(self):
        if not self._app_request("GET"):
            super().do_GET()

    def do_POST(self):
        if not self._app_request("POST"):
            self.send_error(405)

    def translate_path(self, path):
        # mesh manifests "<id>:0" are stored as "<id>_0" (see write_meshes)
        path = re.sub(r":|%3[Aa]", "_", path.split("?", 1)[0].split("#", 1)[0])
        for prefix, folder in self.mounts.items():
            if path.startswith(prefix + "/"):
                saved, self.directory = self.directory, str(folder)
                try:
                    return super().translate_path(path[len(prefix):])
                finally:
                    self.directory = saved
        return super().translate_path(path)

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()

    def log_message(self, *args):
        pass


class QuietHTTPServer(http.server.ThreadingHTTPServer):
    daemon_threads = True
    # Neuroglancer asks for many chunks at once, each on a new connection. With the
    # default listen backlog (5), Windows refuses the extra connections (WinError
    # 10061) instead of letting the client retry like Linux does: those chunks never
    # load and show as empty 64-voxel squares in the slice views.
    request_queue_size = 128

    def handle_error(self, request, client_address):
        # Neuroglancer cancels chunk requests it no longer needs (pan / zoom); the
        # browser then closes the socket mid-transfer. That is normal, not an error.
        if isinstance(sys.exc_info()[1], (BrokenPipeError, ConnectionResetError)):
            return
        super().handle_error(request, client_address)


def open_viewer(out: Path, meta: dict, port: int, layout: str, browser: bool,
                zooms: list | None = None):
    zooms = zooms or []
    mounts = {f"/zoom/{zm['tag']}": zm["dir"] for _, zm in zooms}
    app_ref = {}  # filled below, once the viewer exists
    handler = functools.partial(CORSHandler, directory=str(out), mounts=mounts)

    def make_handler(*a, **kw):
        return handler(*a, app=app_ref.get("app"), **kw)

    server = QuietHTTPServer(("127.0.0.1", port), make_handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    src = f"precomputed://{base}"
    colors = {int(k): v for k, v in meta["colors"].items()}
    initial_luts = {}

    viewer = neuroglancer.Viewer()
    with viewer.txn() as s:
        # Explicit global coordinate space (the registration grid). Without it, zoom
        # values sent from Python (cross_section_scale) are read as metres per pixel.
        res_nm = json.loads((out / "boundaries" / "info").read_text())["scales"][0]["resolution"]
        s.dimensions = neuroglancer.CoordinateSpace(
            names=["x", "y", "z"], units="m", scales=[r * 1e-9 for r in res_nm])
        for i, ch in enumerate(meta["channels"]):
            lut = _LUT_FOR_COLOR[CHANNEL_COLORS[i % len(CHANNEL_COLORS)]]
            initial_luts[ch] = lut
            s.layers[ch] = neuroglancer.ImageLayer(
                source=f"{src}/{ch}", shader=lut_shader(lut),
                shader_controls={"normalized": {"range": meta["ranges"][ch]},
                                 "color": LUTS[lut][2]},
                blend="default" if i == 0 else "additive")
            s.layers[ch].visible = i == 0
        if "template" in meta["images"]:
            s.layers["template"] = neuroglancer.ImageLayer(
                source=f"{src}/template",
                shader_controls={"normalized": {"range": meta["ranges"]["template"]}})
            s.layers["template"].visible = False
        # zoom channels at native resolution, under the atlas outlines / labels
        for zi, (zoom, zm) in enumerate(zooms):
            space = neuroglancer.CoordinateSpace(names=["x", "y", "z"], units="m",
                                                 scales=zoom.atlas_scale_m)
            for ci, ch in enumerate(zm["channels"]):
                color = ZOOM_COLORS.get(ch, CHANNEL_COLORS[(ci + 1) % len(CHANNEL_COLORS)])
                lut = _LUT_FOR_COLOR.get(color, "grays")
                initial_luts[f"{zoom.tag} {ch}"] = lut
                s.layers[f"{zoom.tag} {ch}"] = neuroglancer.ImageLayer(
                    source=neuroglancer.LayerDataSource(
                        url=f"{src}/zoom/{zoom.tag}/{ch}",
                        transform=neuroglancer.CoordinateSpaceTransform(
                            output_dimensions=space, matrix=zoom.matrix)),
                    shader=lut_shader(lut),
                    shader_controls={"normalized": {"range": zm["ranges"][ch]},
                                     "color": LUTS[lut][2]},
                    blend="default" if ci == 0 else "additive")
        if zooms:
            s.position = zooms[0][0].center
        s.layers["boundaries"] = neuroglancer.ImageLayer(
            source=f"{src}/boundaries",
            shader="void main() { float v = float(toRaw(getDataValue()));"
                   " emitRGBA(vec4(1.0, 1.0, 1.0, v > 0.0 ? 0.7 : 0.0)); }")
        for k, lv in enumerate(meta["levels"]):
            s.layers[lv["name"]] = neuroglancer.SegmentationLayer(
                source=f"{src}/{lv['dir']}", selected_alpha=0.35, not_selected_alpha=0.0,
                segment_colors=colors)
            s.layers[lv["name"]].visible = k == 0
        for ex in meta["extra"]:
            s.layers[ex["name"]] = neuroglancer.SegmentationLayer(
                source=f"{src}/{ex['dir']}", selected_alpha=0.25)
            s.layers[ex["name"]].visible = False
        if meta["meshes"]:
            s.layers["region meshes"] = neuroglancer.SegmentationLayer(
                source=f"{src}/meshes", segments=meta["mesh_default"], object_alpha=0.4,
                segment_colors=colors)
        if meta["levels"]:
            first = meta["levels"][0]["name"]
            s.layers[first].tab = "segments"
            s.selected_layer.layer = first
            s.selected_layer.visible = True
        s.layout = layout
        if zooms and layout == "4panel":
            # 4th panel (where the 3D view usually is): an xy view zoomed on the first
            # ROI, showing it at native resolution with the atlas outlines / labels,
            # without the region meshes. Position is shared, zoom level is its own.
            all_layers = [layer.name for layer in s.layers]
            roi_layers = [n for n in all_layers if n != "region meshes"]
            ex, ey, _ = zooms[0][0].extent
            fit = max(ex, ey) * 1.1 / 700  # viewer voxels per screen pixel (~700 px panel)

            def panel(kind, names, **kw):
                return neuroglancer.LayerGroupViewer(layers=names, layout=kind, **kw)

            roi_panel = panel("xy", roi_layers, cross_section_scale=neuroglancer.LinkedZoomFactor(
                link="unlinked", value=fit))
            s.layout = neuroglancer.row_layout([
                neuroglancer.column_layout([panel("xy", all_layers), roi_panel]),
                neuroglancer.column_layout([panel("xz", all_layers), panel("yz", all_layers)]),
            ])

    # Selection panel (right) open: it lists every statistic of the region under the
    # mouse. Double-click pins it there, for that region at every hierarchy level.
    s_json = viewer.state.to_json()
    s_json["selection"] = {}
    viewer.set_state(s_json)
    level_names = [lv["name"] for lv in meta["levels"]]
    seg_names = level_names + [ex["name"] for ex in meta["extra"]]
    if meta["meshes"]:
        seg_names.append("region meshes")  # hit in the 3D view
    app = ViewerApp(viewer, meta, list(initial_luts), level_names, seg_names)
    app.luts = {n: {"lut": lut, "invert": False} for n, lut in initial_luts.items()}
    app.index = RegionIndex(out, meta)

    def select_and_pin(action):
        hits = {}
        for name in seg_names:
            v = action.selected_values.get(name)
            seg = v.value if v is not None else None
            if isinstance(seg, neuroglancer.SegmentIdMapEntry):
                seg = seg.key
            if seg:
                hits[name] = int(seg)
        with viewer.txn() as s:
            hits = {n: v for n, v in hits.items() if s.layers[n].visible}
            if not hits or action.mouse_voxel_coordinates is None:
                s._json_data["selection"] = {}  # double-click on background: unpin
                return
            # usual double-click: show only the clicked region (again: show it no more)
            for name, seg in hits.items():
                segs = s.layers[name].segments
                added = seg not in segs
                if added:
                    segs.add(seg)
                else:
                    segs.discard(seg)
                if name in level_names and "region meshes" in s.layers:
                    msegs = s.layers["region meshes"].segments
                    if added:
                        msegs.add(seg)
                    else:
                        msegs.discard(seg)
            pinned = dict(hits)
            for name, seg in hits.items():
                for k, v in app.pins_for(name, seg).items():
                    pinned.setdefault(k, v)
            s._json_data["selection"] = app._selection(pinned, action.mouse_voxel_coordinates)

    viewer.actions.add("lightsuite-select-and-pin", select_and_pin)
    with viewer.config_state.txn() as c:
        c.input_event_bindings.data_view["dblclick0"] = "lightsuite-select-and-pin"

    ng_url = viewer.get_viewer_url()
    app.viewer_url = ng_url
    app_ref["app"] = app
    url = f"{base}/lightsuite/"
    print(f"\nSample: {meta['sample']}  ({meta['results']})")
    print("Viewer (with region search and channel LUTs):", url)
    print("Plain Neuroglancer:", ng_url)
    if not meta.get("centers"):
        print("  (region centres missing: search highlights regions but cannot jump to them;"
              " run with --reexport)")
    print("Press Ctrl+C to stop.")
    if browser:
        webbrowser.open(url)
    stop = threading.Event()
    try:
        while not stop.wait(0.5):  # short waits so Ctrl+C also works on Windows
            pass
    except KeyboardInterrupt:
        print("Stopping.")
        server.shutdown()


def _writable(d: Path) -> bool:
    """Real write test (os.access is unreliable for folders on Windows / network shares)."""
    if (d / "neuroglancer").is_dir():
        d = d / "neuroglancer"
    probe = d / f".lightsuite_write_test_{os.getpid()}"
    try:
        probe.write_bytes(b"")
        probe.unlink()
        return True
    except OSError:
        return False


def cache_root() -> Path:
    if sys.platform == "win32" and os.environ.get("LOCALAPPDATA"):
        return Path(os.environ["LOCALAPPDATA"]) / "lightsuite-neuroglancer"
    xdg = os.environ.get("XDG_CACHE_HOME")
    return (Path(xdg) if xdg else Path.home() / ".cache") / "lightsuite-neuroglancer"


def output_dir(sample: Sample, out_arg: str | None) -> Path:
    if out_arg:
        return Path(out_arg).expanduser()
    out = sample.root / "neuroglancer"
    if _writable(sample.root):
        return out
    tag = hashlib.sha1(str(sample.root).encode()).hexdigest()[:8]
    return cache_root() / f"{sample.name}-{tag}"


def main():
    for stream in (sys.stdout, sys.stderr):
        try:  # live progress when piped, and no crash on non-UTF-8 Windows consoles
            stream.reconfigure(line_buffering=True, errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(prog="lightsuite-neuroglancer",
                                 description=__doc__.split("\n")[0])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    ap.add_argument("folder", nargs="?", help="LightSuite registration results folder")
    ap.add_argument("--atlas-dir", help="folder with the atlas ontology CSV "
                                        "(default: from the config yaml / ~/atlases)")
    ap.add_argument("--out", help="where to write the converted data "
                                  "(default: <results>/neuroglancer)")
    ap.add_argument("--reexport", action="store_true", help="force re-conversion")
    ap.add_argument("--no-meshes", action="store_true", help="skip 3D region meshes")
    ap.add_argument("--zoom", action="append", default=[], metavar="MULTIRES_DIR",
                    help="LightSuite multires results folder (multires_regopts.json): show its "
                         "raw ROI stacks at native resolution, placed with the saved "
                         "ROI-to-overview transforms (repeatable)")
    ap.add_argument("--port", type=int, default=0, help="data server port (default: any free)")
    ap.add_argument("--layout", default="4panel", help="4panel, xy, yz, xz or 3d (4panel with --zoom: the 3D panel is replaced "
                         "by an xy view zoomed on the ROI, with the atlas labels)")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()

    sample = Sample(choose_results_dir(args.folder))
    out = output_dir(sample, args.out)
    sig = source_signature(sample, args.atlas_dir) + ("|nomesh" if args.no_meshes else "")
    meta_path = out / "viewer_meta.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8")) if meta_path.is_file() else {}
    if args.reexport or meta.get("signature") != sig:
        meta = export(sample, out, args.atlas_dir, meshes=not args.no_meshes)
        meta["signature"] = sig
        meta_path.write_text(json.dumps(meta))
    else:
        print(f"Using cached conversion in {out}")
    zooms = []
    if args.zoom:
        atlas_info = json.loads((out / meta["levels"][0]["dir"] / "info").read_text()) \
            if meta.get("levels") else json.loads((out / "boundaries" / "info").read_text())
        for z in args.zoom:
            zoom = Zoom(Path(z), sample, atlas_info)
            if any(zoom.tag == other.tag for other, _ in zooms):
                zoom.tag = f"{zoom.tag}_{len(zooms) + 1}"
            zooms.append((zoom, prepare_zoom(zoom, args.out, args.reexport)))
            c = ", ".join(f"{v:.1f}" for v in zoom.center)
            print(f"  zoom '{zoom.tag}' ({', '.join(zoom.channels)}) centred at viewer voxel ({c})")
    open_viewer(out, meta, args.port, args.layout, not args.no_browser, zooms)


if __name__ == "__main__":
    main()
