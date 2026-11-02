# /// script
# requires-python = ">=3.10"
# dependencies = ["neuroglancer", "tifffile", "imagecodecs", "numpy", "pandas",
#                 "scipy", "scikit-image", "pyyaml"]
# ///
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

Region names / hierarchy come from, in order of preference:
  1. Fiederling spinal cord atlas Atlas_Regions.csv (parent tree)
  2. Allen parcellation_to_parcellation_term_membership.csv
  3. a downloaded BrainGlobe atlas' structures.json (~/.brainglobe,
     /home/*/.brainglobe or $BRAINGLOBE_DIR), e.g. Perens LSFM, Princeton
  4. the sample's stats/*region_stats.csv (region / structure / division)
  5. bare label numbers
"""

import argparse
import functools
import hashlib
import http.server
import json
import os
import re
import shutil
import struct
import sys
import threading
import webbrowser
from pathlib import Path

import neuroglancer
import numpy as np
import pandas as pd
import tifffile
import yaml
from scipy import ndimage
from skimage.measure import marching_cubes

TOOL_VERSION = 3
MANIFEST = Path("volume_registered") / "sample_space" / "sample_space_manifest.json"
CHUNK = 64
SYNTHETIC_ID_BASE = 3_000_000_000
CHANNEL_COLORS = ["#ffffff", "#00ff00", "#ff00ff", "#00ffff", "#ffff00"]
ALLEN_MEMBERSHIP = "parcellation_to_parcellation_term_membership.csv"
CORD_REGIONS = "Atlas_Regions.csv"


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
        if override:
            dirs.append(Path(override).expanduser())
        cfg_dir = (self.config.get("atlas") or {}).get("atlas_dir")
        if cfg_dir:
            dirs.append(Path(cfg_dir).expanduser())
        if self.manifest.get("segments_csv"):
            dirs.append(Path(self.manifest["segments_csv"]).parent)
        atlases = Path.home() / "atlases"
        if atlases.is_dir():
            dirs += sorted(d for d in atlases.iterdir() if d.is_dir())
        return [d for d in dict.fromkeys(dirs) if d.is_dir()]

    def source_files(self) -> list[Path]:
        files = [self.root / MANIFEST]
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


def fiederling_records(csv: Path) -> list[dict]:
    recs = []
    for _, r in pd.read_csv(csv).iterrows():
        rgb = (int(r["red"]), int(r["green"]), int(r["blue"]))
        recs.append({"id": int(r["id"]), "acronym": str(r["acronym"]),
                     "name": str(r["name"]).strip(),
                     "color": "#%02x%02x%02x" % rgb if any(rgb) else None,
                     "parent": int(r["parent_ID"]) if pd.notna(r["parent_ID"]) else 0})
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
    roots += [Path.home() / ".brainglobe"] + sorted(Path("/home").glob("*/.brainglobe"))
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
                onto = tree_ontology(fiederling_records(d / CORD_REGIONS), labels,
                                     f"atlas region tree ({d / CORD_REGIONS})")
                onto.dir = d
                break
    if onto is None:
        prefer = [d for d in dirs if sample.atlas_name.lower() in d.name.lower()] + dirs
        for d in dict.fromkeys(prefer):
            if (d / ALLEN_MEMBERSHIP).is_file() and sample.atlas_name in ("allen", "unknown"):
                onto = allen_ontology(d / ALLEN_MEMBERSHIP, labels)
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
    return meta


# --------------------------------------------------------------------------
# Viewer
# --------------------------------------------------------------------------

class CORSHandler(http.server.SimpleHTTPRequestHandler):
    def translate_path(self, path):
        # mesh manifests "<id>:0" are stored as "<id>_0" (see write_meshes)
        path = re.sub(r":|%3[Aa]", "_", path.split("?", 1)[0].split("#", 1)[0])
        return super().translate_path(path)

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()

    def log_message(self, *args):
        pass


def open_viewer(out: Path, meta: dict, port: int, layout: str, browser: bool):
    handler = functools.partial(CORSHandler, directory=str(out))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", port), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    src = f"precomputed://http://127.0.0.1:{server.server_address[1]}"
    colors = {int(k): v for k, v in meta["colors"].items()}

    viewer = neuroglancer.Viewer()
    with viewer.txn() as s:
        for i, ch in enumerate(meta["channels"]):
            color = CHANNEL_COLORS[i % len(CHANNEL_COLORS)]
            s.layers[ch] = neuroglancer.ImageLayer(
                source=f"{src}/{ch}",
                shader=f'#uicontrol vec3 color color(default="{color}")\n'
                       "#uicontrol invlerp normalized\n"
                       "void main() { emitRGB(color * normalized()); }",
                shader_controls={"normalized": {"range": meta["ranges"][ch]}})
            s.layers[ch].visible = i == 0
        if "template" in meta["images"]:
            s.layers["template"] = neuroglancer.ImageLayer(
                source=f"{src}/template",
                shader_controls={"normalized": {"range": meta["ranges"]["template"]}})
            s.layers["template"].visible = False
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

    # Selection panel (right) open: it lists every statistic of the region under the
    # mouse. Double-click pins it there, for that region at every hierarchy level.
    s_json = viewer.state.to_json()
    s_json["selection"] = {}
    viewer.set_state(s_json)
    level_names = [lv["name"] for lv in meta["levels"]]
    seg_names = level_names + [ex["name"] for ex in meta["extra"]]
    if meta["meshes"]:
        seg_names.append("region meshes")  # hit in the 3D view

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
                if name in level_names:
                    levels = [level_names.index(name)]
                elif name == "region meshes":
                    levels = range(len(level_names))
                else:
                    continue
                for k in levels:
                    row = next((r for r in meta["hierarchy"] if r[k] == seg), None)
                    if row is not None:
                        for j in range(k, len(level_names)):
                            if row[j]:
                                pinned.setdefault(level_names[j], row[j])
                        break
            order = {n: i for i, n in enumerate(seg_names)}
            s._json_data["selection"] = {
                "position": [float(x) for x in action.mouse_voxel_coordinates],
                "layers": {n: {"value": str(pinned[n])}
                           for n in sorted(pinned, key=order.get)},
            }

    viewer.actions.add("lightsuite-select-and-pin", select_and_pin)
    with viewer.config_state.txn() as c:
        c.input_event_bindings.data_view["dblclick0"] = "lightsuite-select-and-pin"

    url = viewer.get_viewer_url()
    print(f"\nSample: {meta['sample']}  ({meta['results']})")
    print("Neuroglancer viewer:", url)
    print("Press Ctrl+C to stop.")
    if browser:
        webbrowser.open(url)
    try:
        threading.Event().wait()
    except KeyboardInterrupt:
        server.shutdown()


def output_dir(sample: Sample, out_arg: str | None) -> Path:
    if out_arg:
        return Path(out_arg).expanduser()
    out = sample.root / "neuroglancer"
    if os.access(sample.root, os.W_OK):
        return out
    tag = hashlib.sha1(str(sample.root).encode()).hexdigest()[:8]
    return Path.home() / ".cache" / "lightsuite-neuroglancer" / f"{sample.name}-{tag}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder", nargs="?", help="LightSuite registration results folder")
    ap.add_argument("--atlas-dir", help="folder with the atlas ontology CSV "
                                        "(default: from the config yaml / ~/atlases)")
    ap.add_argument("--out", help="where to write the converted data "
                                  "(default: <results>/neuroglancer)")
    ap.add_argument("--reexport", action="store_true", help="force re-conversion")
    ap.add_argument("--no-meshes", action="store_true", help="skip 3D region meshes")
    ap.add_argument("--port", type=int, default=0, help="data server port (default: any free)")
    ap.add_argument("--layout", default="4panel", help="4panel, xy, yz, xz or 3d")
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
    open_viewer(out, meta, args.port, args.layout, not args.no_browser)


if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    main()
