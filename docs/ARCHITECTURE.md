# Architecture and maintenance guide

For developers who maintain or extend `lightsuite-neuroglancer`. For installing and
using the tool, see [installation.html](installation.html) and
[user-guide.html](user-guide.html).

Contents:

1. [What the tool does](#1-what-the-tool-does)
2. [Repository layout](#2-repository-layout)
3. [Dependencies](#3-dependencies)
4. [Runtime overview](#4-runtime-overview)
5. [Input contract (what LightSuite must produce)](#5-input-contract-what-lightsuite-must-produce)
6. [Conversion pipeline](#6-conversion-pipeline)
7. [Output layout and cache](#7-output-layout-and-cache)
8. [Coordinates and orientation](#8-coordinates-and-orientation)
9. [Zoom stacks](#9-zoom-stacks)
10. [Viewer, side panel and HTTP API](#10-viewer-side-panel-and-http-api)
11. [Packaging and releases](#11-packaging-and-releases)
12. [Maintenance notes and known issues](#12-maintenance-notes-and-known-issues)
13. [How to extend](#13-how-to-extend)
14. [Manual test checklist](#14-manual-test-checklist)

---

## 1. What the tool does

It reads the output folder of a LightSuite atlas registration and shows it in
[Neuroglancer](https://github.com/google/neuroglancer). It works in two phases:

1. **Convert once.** The TIFF volumes, atlas ontology, region statistics and
   meshes are turned into Neuroglancer's
   [precomputed format](https://github.com/google/neuroglancer/tree/master/src/datasource/precomputed)
   in `<results>/neuroglancer/`. The result is cached and reused.
2. **Serve and view.** A local HTTP server serves the converted files. A
   Neuroglancer Python viewer is configured with the layers, and a small side
   panel (region search, channel LUTs) is wrapped around it in the browser.

Design choices worth keeping:

- **One file.** All the code is in `view_in_neuroglancer.py`, with an inline
  [PEP 723](https://peps.python.org/pep-0723/) dependency header, so
  `uv run view_in_neuroglancer.py` works with no install step. The same file is also
  the installable module (`py-modules` in `pyproject.toml`).
- **Read-only on the inputs.** The tool never changes LightSuite files. It writes only
  to `neuroglancer/` (or to a cache folder, see §7).
- **No re-registration, no resampling.** Everything is shown on the pipeline's own
  registration grid. Zoom stacks are placed with one affine built from transforms the
  pipeline saved earlier.
- **Local only.** All servers bind to `127.0.0.1`.

## 2. Repository layout

| Path | Role |
|---|---|
| `view_in_neuroglancer.py` | The whole application (about 2,000 lines). Entry point `main()`. |
| `pyproject.toml` | Packaging metadata, dependencies, console scripts. |
| `run_viewer.bat` | Windows launcher: checks that `uv` exists, then `uv run`s the script with the dropped folder or arguments. |
| `README.md` | Short overview, options, data sources. |
| `docs/installation.html` | End-user installation guide. |
| `view_in_neuroglancer.py.lock` | Exact, tested dependency versions for `uv run view_in_neuroglancer.py` (uv script lock file). Commit it. |
| `docs/user-guide.html` | End-user usage guide (tasks and interpretation). `docs/images/` holds its screenshot. |
| `docs/troubleshooting.html` | End-user troubleshooting: run checklist, console output, every error message. |
| `docs/ARCHITECTURE.md` | This file. |
| `.gitignore` | Ignores build output, `.venv/` and the project `uv.lock` (the script lock above is not ignored). |

There are no tests, CI or subpackages.

### Sections of `view_in_neuroglancer.py`

The file is split into sections with `# ----` banner comments, in this order:

| Section | Main symbols | Purpose |
|---|---|---|
| Header | PEP 723 block, module docstring, constants | `__version__`, `TOOL_VERSION`, `MANIFEST`, `CHUNK`, `SYNTHETIC_ID_BASE`, colours. |
| Locating results | `find_results_dirs`, `ask_folder`, `choose_results_dir` | Find the folder(s) with `volume_registered/sample_space/sample_space_manifest.json` (up to 4 levels deep). Folder picker with tkinter, or a terminal prompt as fallback. |
| Sample metadata | `Sample` | Loads the manifest, config YAML, stats CSVs. Resolves paths that moved (`Sample.resolve`). Lists atlas folders and source files. |
| Volumes | `load_yxz`, `permute_brain_volume`, `to_xyz`, `compute_boundaries`, `load_volumes` | Reads TIFFs and puts them in viewer orientation. |
| Ontology | `Ontology`, `allen_ontology`, `tree_ontology`, `stats_ontology`, `build_ontology`, … | Maps each label to one term per hierarchy level. |
| Per-region statistics | `attach_metrics`, `ID_COLUMNS` | Adds numeric properties from LightSuite CSVs. |
| Precomputed writers | `write_precomputed`, `write_segment_properties`, `write_meshes`, `region_centers`, `extra_names` | Writes the Neuroglancer data formats. |
| Export | `source_signature`, `export` | Runs the whole conversion and returns `meta`. |
| Zoom | `Zoom`, `_elastix_*`, `_registration_grid_matrix`, `write_zoom_channel`, `prepare_zoom` | Multires ROI stacks (§9). |
| LUTs | `_TABLES`, `LUTS`, `lut_rgb`, `lut_css`, `lut_shader` | Colour maps as GLSL shaders and CSS gradients. |
| Region search | `_norm`, `RegionIndex` | Fuzzy search index over all segment layers. |
| Control panel | `ViewerApp`, `PANEL_HTML` | Server-side panel actions and the panel page (HTML/CSS/JS kept in one string). |
| HTTP | `CORSHandler`, `QuietHTTPServer` | Static file server, plus the `/lightsuite/` panel and API. |
| Viewer | `open_viewer` | Builds the Neuroglancer state and runs until Ctrl+C. |
| Output location | `_writable`, `cache_root`, `output_dir` | Where converted data goes. |
| CLI | `main` | Argument parsing and orchestration. |

## 3. Dependencies

### Python packages (runtime)

Declared twice, **keep them in sync**: in the `# /// script` header of
`view_in_neuroglancer.py` and in `[project].dependencies` of `pyproject.toml`.
Python ≥ 3.10 is required (the code uses `X | None`, `str.removeprefix`, the
walrus operator and the `dict | dict` merge).

Versions are fixed at two levels:

| Level | Where | Used by |
|---|---|---|
| **Ranges** | both dependency lists: `neuroglancer==2.41.*` (one minor version, see §12), `numpy>=1.24,<3`, `pandas>=2.0,<4`, `scipy>=1.10,<2`, `scikit-image>=0.21,<1`, `pyyaml>=6,<7`; `tifffile` and `imagecodecs` (calendar versions) have lower bounds only | `uv tool install`, `pip install` |
| **Exact lock** | `view_in_neuroglancer.py.lock`, made by `uv lock --script view_in_neuroglancer.py` | `uv run view_in_neuroglancer.py` and `run_viewer.bat` |

The lock resolves for every Python ≥ 3.10 at once, so it holds the newest versions
that still support 3.10 (in November 2026: numpy 2.4, scipy 1.17, pandas 3.1,
scikit-image 0.26, neuroglancer 2.41.2). It was tested with a full conversion and
viewer run on Python 3.11 (Windows). The lower bounds of the ranges were not tested.

| Package | Used for | Where |
|---|---|---|
| `neuroglancer` | Python viewer: runs its own web server (Tornado) that serves the Neuroglancer JS client and syncs the viewer state over a websocket. Also provides the state classes (`ImageLayer`, `SegmentationLayer`, `CoordinateSpace`, layouts, actions). | `open_viewer`, `ViewerApp` |
| `numpy` | All volume handling, LUT maths, affine matrices. | everywhere |
| `tifffile` | Reads multi-page TIFFs (`TiffFile`, `memmap` for the large zoom stacks). | `load_yxz`, `_load_stack_zyx` |
| `imagecodecs` | Not imported directly. `tifffile` needs it to decode compressed TIFFs (LZW, Deflate, …). | indirect |
| `pandas` | Reads the ontology and statistics CSVs. | ontology, `attach_metrics` |
| `scipy` (`ndimage`) | `find_objects` (region bounding boxes), `gaussian_filter` (smoothing before meshing), `center_of_mass`. | `write_meshes`, `region_centers` |
| `scikit-image` | `marching_cubes` for region surface meshes. | `write_meshes` |
| `pyyaml` | Reads the LightSuite sample config (`*.yaml` / `*.yml` in the results folder). | `Sample._read_config` |

### Standard library worth knowing about

`http.server` (data and panel server), `threading` (server thread, panel lock),
`hashlib` (cache signatures), `secrets` (panel API token), `difflib` and
`unicodedata` (fuzzy region search), `webbrowser`, `tkinter` (optional folder
picker; on Linux without `python3-tk` or a display the tool falls back to a terminal
prompt), `concurrent.futures` (parallel zoom chunk writing).

### Outside Python

| Item | Why |
|---|---|
| A browser with WebGL 2 | The Neuroglancer client renders with WebGL 2. |
| `uv` (recommended) | Gets Python and the dependencies from the PEP 723 header or `pyproject.toml`. `pip` also works. |
| Internet on the first run | Only to download the packages. The Neuroglancer client comes inside the `neuroglancer` wheel, so viewing works offline. |
| Atlas files (optional) | Allen membership CSV, BrainGlobe `structures.json`, or the Fiederling `Atlas_Regions.csv`, for region names (§6.3). |

## 4. Runtime overview

```
                 ┌────────────────────────── python process ───────────────────────────┐
                 │                                                                      │
 CLI / .bat ──►  main()                                                                 │
                 │  choose_results_dir ─► Sample ─► output_dir                          │
                 │        │                                                             │
                 │        ▼  signature changed?                                         │
                 │  export()  ──writes──►  <results>/neuroglancer/ + viewer_meta.json   │
                 │        │                                                             │
                 │  Zoom / prepare_zoom (per --zoom) ──► <multires>/neuroglancer/       │
                 │        │                                                             │
                 │  open_viewer()                                                       │
                 │    ├── QuietHTTPServer (port --port, thread)                         │
                 │    │      GET  /<layer>/...         precomputed files (CORS *)       │
                 │    │      GET  /zoom/<tag>/...      zoom precomputed files           │
                 │    │      GET  /lightsuite/         side panel page                  │
                 │    │      GET|POST /lightsuite/api/* panel API (token)               │
                 │    └── neuroglancer.Viewer (own Tornado server, other port)          │
                 │           /v/<id>/  Neuroglancer client + websocket state            │
                 └──────────────────────────────────────────────────────────────────────┘

 Browser:  http://127.0.0.1:<port>/lightsuite/
           ┌─ side panel (PANEL_HTML) ─┐┌─ <iframe> Neuroglancer client (/v/<id>/) ─┐
           │ search, LUTs ──fetch──► API ││  loads precomputed:// data from <port>    │
           └────────────────────────────┘└───────────────────────────────────────────┘
```

Two servers run at the same time:

- **Data and panel server** (`QuietHTTPServer` + `CORSHandler`), on `--port` (by
  default a free port). It serves the converted folder, mounts zoom folders
  under `/zoom/<tag>`, and handles `/lightsuite/…` through `ViewerApp.handle`.
  It sends `Access-Control-Allow-Origin: *` because the Neuroglancer client is
  served from the other port.
- **Neuroglancer's own server** (started by `neuroglancer.Viewer()`). It serves the
  client and keeps the Python-side `viewer.state` in sync with the browser.
  Panel actions change `viewer.state` in Python (`viewer.txn()`), and the change is
  sent to the browser over this websocket.

The main thread waits in a `stop.wait(0.5)` loop so that Ctrl+C is handled on
Windows too. Then it shuts the data server down. Viewer state is **not** kept
between runs.

## 5. Input contract (what LightSuite must produce)

The tool depends on these LightSuite outputs. Changes on the LightSuite side need to
be checked against this list.

### Required

`volume_registered/sample_space/sample_space_manifest.json`, with:

| Key | Meaning |
|---|---|
| `annotation_path` | Label volume (TIFF) on the registration grid. **Required.** |
| `channel_paths` | `{ "1": path, "2": path, … }` sample channels. |
| `voxel_um` | Registration-grid voxel size, in `(y, x, z)` order (default `[20, 20, 20]`). |
| `permute_sample_to_atlas` | Axis permutation/flip vector (1-based, negative = flip). Default `[1, 2, 3]`. |
| `atlas_volumes_permuted` | `true` if atlas volumes were already written permuted. |
| `template_path`, `boundary_path` | Optional. Boundaries are computed from the annotation if missing. |
| any other `*_path` | Shown as an extra label layer (e.g. `segments_path`, `division_labels_path`) if it is not empty. |
| `segments_csv` | Spinal cord: segment names; its folder is also searched for `Atlas_Regions.csv`. |

All volumes must have the annotation's shape, otherwise the tool stops.
Paths in the manifest may be absolute paths from another machine. `Sample.resolve`
then looks for the same file name in `sample_space/` and in the results root.

### Optional, used when present

| File | Used for |
|---|---|
| `<results>/*.yaml` / `*.yml` (first one with an `atlas:` mapping) | `atlas.provider`, `atlas.atlas_dir`, `atlas.brainglobe_name`, `atlas.resolution_um`, `sample.name` |
| `regopts.json` | `brain_atlas` (atlas name fallback). For `--zoom`: `nx, ny, nz, voxel_um, registres_um`. |
| `registration_analysis/region_metrics.csv` | Statistics at every hierarchy level (columns `channel`, `level`, `region_id`, `acronym` + numeric columns). |
| `stats/region_stats.csv` or `stats/chanNN_region_stats.csv` | Atlas-space statistics (long format: `parcellation_index`, `channel`, `hemisphere`, `metric`, `value`). Also a names fallback. |
| `stats/region_stats_sample.csv` or `stats/chanNN_region_stats_sample.csv` | Sample-space statistics, same format. |
| `division_id_legend*.csv` in the atlas folder | Names for the `division_labels` extra layer. |

## 6. Conversion pipeline

`export()` runs these steps in order. It always starts from an empty output folder
(`shutil.rmtree(out)`).

### 6.1 Volumes (`load_volumes`)

1. `load_yxz` stacks TIFF pages into a `(Y, X, Z)` array, as LightSuite's
   `load_registration_volume` does.
2. Channels always get `permute_brain_volume(v, perm)`. Atlas volumes (annotation,
   template, boundaries, extras) get it only if `atlas_volumes_permuted` is false.
3. `to_xyz` swaps the first two axes, so the arrays are `(x, y, z)`.
4. The voxel size in nm is `[vox[1], vox[0], vox[2]] * 1000`.

### 6.2 Image layers

`write_precomputed` writes one scale (`s0`), `raw` encoding, 64³ chunks, in
Fortran order (x varies fastest). Segmentations are `uint32`; images keep their
dtype (booleans become `uint8`). The default contrast range is the 1st–99.5th
percentile of the non-zero voxels (every 2nd voxel per axis), saved in
`meta["ranges"]`.

### 6.3 Ontology (`build_ontology`)

The `Ontology` object holds:

- `levels`: hierarchy level names, **finest first**;
- `chains`: `label → [term key per level]` (`None` = no parent at that level);
- `terms`: `term key → {acronym, name, color, native_id, metrics}`;
- `ids`: `term key → segment id` (`assign_ids`). The atlas's own ID is kept when it
  is unique and `< 3,000,000,000`. Other terms get synthetic IDs from
  `SYNTHETIC_ID_BASE`.

Providers are tried in this order. The first one that matches wins:

| # | Provider | Condition | Builder |
|---|---|---|---|
| 1 | Fiederling spinal cord `Atlas_Regions.csv` | manifest has `segments_csv` or atlas is `fiederling` | `tree_ontology(table_records)` |
| 2 | Allen `parcellation_to_parcellation_term_membership.csv` | atlas `allen`/`unknown`, ≥ 90 % of labels known | `allen_ontology` |
| 3 | Gubra / Perens LSFM `ARA2_annotation_info.csv` (`id, parent_id, name, acronym, red, green, blue`) | ≥ 90 % of labels known | `tree_ontology(table_records)` |
| 4 | BrainGlobe `structures.json` | ≥ 90 % of labels known; folders ranked by `brainglobe_name`, atlas name, resolution | `tree_ontology(brainglobe_records)` |
| 5 | Sample `stats/region_stats*.csv` | always, if present | `stats_ontology` (levels `regions/structure/division`) |
| 6 | Label numbers only | fallback | empty `Ontology` |

Atlas folders are searched in this order: `--atlas-dir` and its subfolders, config
`atlas.atlas_dir` and its subfolders, the folder of `segments_csv`, `~/atlases/*`.
For providers 2 and 3, folders whose name contains the atlas name (config
`atlas.provider`) come first. BrainGlobe folders: `$BRAINGLOBE_DIR`, `~/.brainglobe`,
other users' `*/.brainglobe` (next to the home folder and in `/home`).

For tree atlases (1, 3, 4) the levels are **depths counted from the root** (at most
4 coarse levels plus the leaf labels), so that each layer cuts the tissue at one
consistent depth. Parents with no colour get the mean colour of their leaves.
`complete()` adds a `label N (not in atlas)` term for labels the provider does not
know.

`warn_label_overflow` then prints a warning when the atlas has IDs larger than the
annotation's integer type can hold and some labels are unknown. This catches
annotations that were saved as `uint16` (often after a float32 step) by the
registration pipeline: large Allen-style IDs (e.g. 312782624) are rounded and
wrapped, and several regions end up sharing one meaningless label. The viewer cannot
repair this; the fix belongs in the pipeline.

### 6.4 Statistics (`attach_metrics`)

- `region_metrics.csv`: every numeric column that is not an ID column (`ID_COLUMNS`)
  and that varies within a channel. Columns that vary between channels get a
  `_chN` suffix. Rows are matched by `level` + `acronym` (or `region_id` at the
  `label` level).
- `region_stats_sample.csv` (no prefix) and `region_stats.csv` (prefix `atlas_`):
  long-format rows become properties named `<prefix><metric>[_chN]_<hemisphere>`.
  Volume metrics never get a channel suffix. Tables with duplicate keys (spinal cord
  per-segment tables) are skipped.
- Missing values are written as the string `"NaN"`, which Neuroglancer hides.
- Each property has a description (shown on hover) that names its source file.

### 6.5 Segmentation layers

- One layer `level_<tag>` per hierarchy level that **actually splits the tissue**
  (more than one term). Made with a lookup table `lut[annot]`, label → term id.
- `write_segment_properties` adds inline properties: `label` (acronym),
  `description` (name), `tags` (the coarser regions it belongs to plus the level tag,
  searchable as `#Isocortex`), and one `number` property per statistic.
- Extra label volumes become `extra_<key>` layers. `division_labels` is skipped
  when a `division` level already exists.
- `meta["hierarchy"]`: one tuple of segment ids (one per level) per atlas label. The
  viewer uses it to pin a region together with all the regions containing it.

### 6.6 Meshes (`write_meshes`) and centres (`region_centers`)

- One mesh per term (every level). The region mask is cropped to its bounding box,
  smoothed with a Gaussian (σ = 1 voxel), and meshed with `marching_cubes` at 0.5
  (step 2 for masks over 8 M voxels). Vertices are in nm.
- Neuroglancer's **legacy mesh format**: manifest `<id>:0` → fragment `<id>_0_1`.
  Windows does not allow `:` in file names, so manifests are saved as `<id>_0` and
  `CORSHandler.translate_path` maps `:` / `%3A` to `_`.
- `mesh_default` = the largest term at the coarsest level (the whole brain or cord),
  shown at start.
- `region_centers`: for each term, a voxel inside the region near its centroid
  (for non-convex regions, the region voxel nearest the centroid) plus the bounding
  box size. The panel uses them to jump to a region and zoom to fit it.

Meshing is the slowest step (single-threaded). `--no-meshes` skips it.

## 7. Output layout and cache

```
<results>/neuroglancer/
├── viewer_meta.json          # meta returned by export() + "signature"
├── ch1/, ch2/, …             # image layers: info + s0/<x0-x1_y0-y1_z0-z1>
├── template/, boundaries/
├── level_<tag>/              # segmentation per hierarchy level
│   ├── info, s0/…
│   └── segment_properties/info
├── extra_<key>/              # extra label volumes
└── meshes/                   # info, <id>_0 manifests, <id>_0_1 fragments, segment_properties/
```

`viewer_meta.json` keys: `version`, `sample`, `results`, `channels`, `images`,
`ranges`, `levels` (`[{name, dir}]`), `extra`, `hierarchy`, `mesh_default`,
`meshes`, `colors` (`segment id → hex`), `centers`, `signature`.

**Where it goes** (`output_dir`): `--out` if given. Otherwise `<results>/neuroglancer`
if a real test write works (`_writable`; `os.access` is not reliable on Windows and
network shares). Otherwise
`cache_root()/<sample>-<sha1(root)[:8]>`, where `cache_root()` is
`%LOCALAPPDATA%\lightsuite-neuroglancer` on Windows, or
`$XDG_CACHE_HOME`/`~/.cache` + `lightsuite-neuroglancer`.

**When it is rebuilt.** `source_signature` is a SHA-1 of:

- `TOOL_VERSION` and `--atlas-dir`;
- **the bytes of `view_in_neuroglancer.py` itself**, so any edit to the file
  (even a comment) causes a full re-conversion on the next run;
- path, size and mtime of every source file (`Sample.source_files`: manifest,
  config YAML, all `*_path` volumes, channels, statistics CSVs);
- `|nomesh` when `--no-meshes` is used.

The atlas ontology files themselves (Allen CSV, Perens CSV, BrainGlobe JSON) are
**not** part of the signature: which one is used is only known during the
conversion. After installing an atlas without changing the config or the
`--atlas-dir` option, the user must run with `--reexport`.

## 8. Coordinates and orientation

- The viewer's global coordinate space is set explicitly in `open_viewer`:
  dimensions `x, y, z`, unit `m`, scales = registration-grid voxel size. Without it,
  `cross_section_scale` values set from Python would be read as metres per pixel.
- Viewer voxel `i` spans `[i, i+1)`. Its centre is `i + 0.5` (this matters for
  `region_centers` and for the zoom affine).
- The orientation is the same as LightSuite's QC overlays: the registration grid
  after `permute_sample_to_atlas`. The functions `load_yxz` and
  `permute_brain_volume` copy LightSuite's `load_registration_volume` and
  `permute_brain_volume` / `permuteBrainVolume.m`. **If LightSuite changes these, the
  copies here must be changed too.**

## 9. Zoom stacks

`--zoom DIR` (repeatable) adds the raw ROI channels of a LightSuite **multires** run.

1. `Zoom.__init__` finds `multires_regopts.json` (in `DIR` or one level below), the
   pair manifest (`converted/*_pair.json`), the ROI TIFF of each channel and the
   elastix `TransformParameters.*.txt` files. All paths go through `_resolve_near`,
   so moved data still works. All channels must have the same ROI geometry.
2. The placement is a single 4×4 affine, ROI voxel → atlas viewer voxel:

   ```
   M = G · inv(I_overview) · inv(E) · L · I_roi
   ```

   | Term | Source |
   |---|---|
   | `I_roi`, `I_overview` | `_index_to_physical(spec)`: `direction · diag(spacing) + origin` from the pair manifest |
   | `L` | `roi_to_overview_tform` (landmark fit) in `multires_regopts.json`; identity if missing (metadata-only geometry) |
   | `E` | `_elastix_chain`: the last TransformParameters file composed with its `InitialTransformParameterFileName` chain. Fixed → moving. Only Translation, Euler and Affine are supported; **B-spline stops the tool**. |
   | `G` | `_registration_grid_matrix`: the downsampling LightSuite's `preprocess.slice_ops` uses for the registration grid (integer stride + corner-aligned zoom in XY, linear in Z), then the `permute_sample_to_atlas` permutation/flip and the x/y swap |

   It is then shifted by ±0.5 for voxel centres. The overview stack in the pair
   manifest must be the stack the atlas run used (`regopts.json` size and voxel
   size), otherwise the tool stops.
3. `write_zoom_channel` writes a multiscale image pyramid: `s0` = native voxels, each
   next level is a 2×2×2 mean (`_halve`), until the largest side is ≤ 128. Chunks
   are written in parallel threads. `info` is written last, so it marks a complete
   conversion.
4. `prepare_zoom` caches in `<multires>/neuroglancer/` (or `--out/zoom_<tag>`, or
   `cache_root()/zoom_<tag>-<hash>` if the multires folder is read-only, with the same
   `_writable` test as the atlas output). The zoom
   signature hashes `ZOOM_VERSION`, `CHUNK`, and each channel's TIFF path/size/mtime/
   spacing. It does **not** hash the script, so bump `ZOOM_VERSION` when you change
   the zoom output format. The affine is recomputed every run and written to
   `zoom_meta.json` for reference.
5. In `open_viewer`, each zoom channel is an `ImageLayer` whose `LayerDataSource`
   carries the affine (`CoordinateSpaceTransform`). With the `4panel` layout the
   3D panel is replaced by an extra xy panel with its own zoom
   (`LinkedZoomFactor(link="unlinked")`) that fits the first ROI.

## 10. Viewer, side panel and HTTP API

### Layer set-up (`open_viewer`)

| Layer | Type | Default |
|---|---|---|
| `ch1…chN` | Image, LUT shader, `ch2+` blended additively | only `ch1` visible |
| `template` | Image | hidden |
| zoom channels `<tag> <channel>` | Image with affine | visible |
| `boundaries` | Image, white at 70 % where > 0 | visible |
| one per hierarchy level | Segmentation, `selected_alpha 0.35`, `not_selected_alpha 0` | only the finest visible |
| extra label layers | Segmentation | hidden |
| `region meshes` | Segmentation (meshes only), `object_alpha 0.4` | `mesh_default` |

The **Selection** panel is opened by setting `selection = {}` in the state JSON. The
mouse binding `dblclick0` is mapped to the Python action
`lightsuite-select-and-pin` (`select_and_pin`): it toggles the region in the layer's
visible segments and in `region meshes`, and pins the region and its parents
(`ViewerApp.pins_for`) in the Selection panel.

### Side panel

`GET /lightsuite/` returns `PANEL_HTML`, with `/*CONFIG*/null` replaced by a JSON
config (API token, Neuroglancer URL, sample name, LUT list with CSS previews,
number of regions). The page shows the Neuroglancer viewer in an `<iframe>`
and calls the API with `fetch`.

| Endpoint | Method | Body / query | Action |
|---|---|---|---|
| `/lightsuite/api/search` | GET | `q` | `RegionIndex.search` (max 25 results) |
| `/lightsuite/api/state` | GET | | channel visibility/LUT and current highlight (the panel polls it every 3 s) |
| `/lightsuite/api/show` | POST | `layer, id, add` | `ViewerApp.show_region`: isolate/highlight, jump, zoom to fit, show mesh, pin |
| `/lightsuite/api/clear` | POST | | restore visibility, alpha, position and zoom saved before the first search |
| `/lightsuite/api/lut` | POST | `layer, lut, invert` | replace the shader, keeping the `normalized` contrast range |
| `/lightsuite/api/visible` | POST | `layer, visible` | toggle a channel |

All API calls need the per-run random `token` (`secrets.compare_digest`).
`ViewerApp.lock` makes state edits from concurrent requests run one at a time.

### Region search (`RegionIndex`)

The index is built from the `segment_properties/info` files of every level and extra
layer. Duplicates (same acronym and name at a coarser level) are dropped, keeping the
finest one. `_norm` folds case, accents, punctuation and `ae`/`oe` → `e`. Scores:
exact 100, name prefix 92, acronym prefix 88, word prefix 84, substring 78, all words
prefixed 74, then fuzzy (`difflib`) up to 70 (shown as "similar" in the panel). A
numeric query matches segment IDs exactly.

### LUTs

`LUTS` maps an id to `(label, kind, parameter)`. Kinds are `color` (black → colour,
with a `#uicontrol vec3 color`), `table` (33 samples of a matplotlib colormap,
interpolated in GLSL), and `fn` (Fire, Jet, HiLo as formulas). `lut_shader` writes the
GLSL. `lut_rgb` runs the same formulas in numpy for the CSS previews. **Keep the two
in sync.**

## 11. Packaging and releases

- `pyproject.toml`: setuptools, `py-modules = ["view_in_neuroglancer"]`, two console
  scripts (`lightsuite-neuroglancer`, `view-in-neuroglancer`) → `main`.
- Three ways to install: `uv run` the file (PEP 723 header), `uv tool install`
  (git URL or local folder), `pip install .`.

### Release checklist

1. Change the version in **five places**: `__version__` in the script,
   `version` in `pyproject.toml`, and the footers of `docs/installation.html`,
   `docs/user-guide.html` and `docs/troubleshooting.html`.
2. If dependencies changed, change **both** the PEP 723 header and `pyproject.toml`,
   then refresh the lock: `uv lock --script view_in_neuroglancer.py` (add
   `--upgrade` to move to the newest allowed versions). Commit the `.lock` file.
3. Bump `TOOL_VERSION` when the **meaning** of the converted output changes (the
   script hash already forces re-conversion on any edit, but the version number
   documents format changes and is saved in `viewer_meta.json`).
4. Bump `ZOOM_VERSION` when the zoom output changes (its cache does not hash the
   script).
5. Run the manual checks in §14, on Windows and on Linux if possible.
6. Update `README.md` and the docs for new options or behaviour.

## 12. Maintenance notes and known issues

**Dependency risk**

- `neuroglancer` is pinned to one minor version (`==2.41.*`) because the code uses
  internal or less stable parts of its API:
  - `s._json_data["selection"]` (private attribute) to open and pin the Selection
    panel;
  - `neuroglancer.SegmentIdMapEntry` in action selected values;
  - `LinkedZoomFactor`, `LayerGroupViewer`, `row_layout`/`column_layout`,
    `config_state.input_event_bindings`, `status_messages`.

  To move to a newer Neuroglancer: change the pin in both lists, `uv lock --script
  --upgrade-package neuroglancer view_in_neuroglancer.py`, and run the whole §14
  checklist (the double-click pin, the panel search and the zoom panel are what break
  first).
- The other packages have upper bounds at their next major version.
  `marching_cubes`, `tifffile.memmap` and the pandas groupby calls are stable APIs, but
  check them before raising a bound.
- `uv tool install` and `pip install` use the ranges, not the lock, so they can get
  newer (untested) minor versions than `uv run`. If a user reports a problem that you
  cannot reproduce, compare versions first.

**Memory and speed**

- The whole registration-grid volume set is loaded into RAM (about 8 GB needed for a
  20 µm brain, including meshing). Zoom stacks are read fully into memory too
  (`np.array(tifffile.memmap(...))`), then the pyramid is built in memory.
- Meshing runs one region at a time in one thread. It is the main cost of the first
  run.
- Atlas layers have a single scale (no pyramid). This is fine at 10–25 µm, but would
  need multiscale output for finer atlas grids.

**Inconsistencies to clean up**

- Atlas ontology files are not part of the cache signature (§7).
- For the Perens table and BrainGlobe atlases, the hierarchy levels are tree depths
  (`depth 3`, `depth 2`, …). These are coarser than LightSuite's own
  `structure`/`division` grouping in `stats/region_stats.csv` (e.g. 39 and 15
  regions instead of 311 and 23 for a Perens brain). Using the stats columns as the
  levels, and the atlas table only for names, colours and parents, would match
  LightSuite's analysis better.

**Upstream data issues the viewer detects**

- Annotation saved as `uint16` with Allen-style IDs (§6.3, `warn_label_overflow`).
  Seen in `8_DNF_atlased` (`annotation_in_sample_20um.tif`): 128 Perens regions folded
  into 23 labels. Belongs in the LightSuite registration export.

**Platform details (do not "simplify" these away)**

- `:` in mesh manifest names (§6.6).
- `stream.reconfigure(line_buffering=True, errors="replace")` in `main`: live progress
  when output is piped, and no crash on non-UTF-8 Windows consoles.
- `stop.wait(0.5)` loop instead of `join()`, so that Ctrl+C works on Windows.
- `QuietHTTPServer.request_queue_size = 128`: with the default backlog of 5,
  Windows refuses connections when Neuroglancer requests many chunks at once
  (`WinError 10061`). Those chunks never load and appear as empty 64-voxel squares in
  the slice views. Linux makes the client retry instead, so the bug is Windows-only.
  A stress test with 64 parallel clients gave 2–9 failures out of 2,178 requests
  before the fix and none after.
- `QuietHTTPServer.handle_error` ignores `BrokenPipeError`/`ConnectionResetError`:
  Neuroglancer cancels chunk requests while you pan.
- Text files are read with `encoding="utf-8"` everywhere. Keep doing that (the Windows
  default code page would break on accented names).

**Not supported yet**

- Spinal cord per-segment statistics tables (skipped in `attach_metrics`).
- B-spline elastix stages for zoom stacks.
- Saving or sharing the viewer state between runs.
- Access from another machine (servers bind to `127.0.0.1` only).

## 13. How to extend

| Task | Where |
|---|---|
| New CLI option | `main()` argparse. Pass it on to `export`/`open_viewer`. If it changes the converted output, add it to the signature in `main` (as `--no-meshes` does). Document it in `README.md` and `docs/installation.html`. |
| New atlas provider | If it is a CSV with `id, parent_id, name, acronym, red, green, blue`, `table_records` already reads it: only add the file name to `build_ontology`. Otherwise write a `*_records` function that returns `[{id, acronym, name, color, parent}]` and use `tree_ontology`, or build an `Ontology` directly. Add it to `build_ontology` in the right priority position. |
| New statistics source | `attach_metrics`: put values in `onto.terms[key]["metrics"][name]` and add `names[name] = description`. Add the file to `Sample.source_files` so the cache is invalidated when it changes. |
| New LUT | Add an entry to `LUTS` (and `_TABLES` for a sampled colormap). For `fn`, add the formula in **both** `lut_rgb` and `lut_shader`. |
| New panel action | Add a method to `ViewerApp`, a branch in `ViewerApp.handle`, and the UI/JS in `PANEL_HTML`. Edit the state inside `with self.lock, self.viewer.txn() as s:`. |
| New layer type | Write it in `export` (add to `meta`), then create the layer in `open_viewer`. If it is a segmentation layer, add it to `seg_names` so hover/pin works. |
| New manifest volume | Any extra `*_path` key is already shown as an extra label layer. Add a naming rule in `extra_names` if needed. |

## 14. Manual test checklist

There are no automated tests. Before a release, run with at least one Allen brain,
one BrainGlobe brain and (if available) one spinal cord dataset:

- [ ] `lightsuite-neuroglancer --version` prints the new version.
- [ ] First run converts. Second run prints `Using cached conversion`.
- [ ] `--reexport` and `--no-meshes` work. Changing a source file triggers re-conversion.
- [ ] Region names and hierarchy levels are correct (the console line `region names from:` shows the expected provider).
- [ ] Hover shows names. The Selection panel shows statistics with descriptions.
- [ ] Double-click isolates, shows the mesh and pins the parents. Double-click again restores. Double-click the background unpins.
- [ ] Panel search: exact name, acronym, ID, a typo, `coeruleus`. Enter jumps and zooms. Shift+click adds. *Show all regions* restores the view.
- [ ] LUT change keeps the contrast range. Invert works. HiLo shows blue/red at the range ends.
- [ ] Orientation matches LightSuite's QC overlays (channel vs boundaries).
- [ ] `--zoom` places the ROI on the matching anatomy. The 4th panel shows the ROI.
- [ ] Read-only results folder: output goes to the cache folder.
- [ ] Windows: `run_viewer.bat` by double-click and by dropping a folder on it. Ctrl+C stops cleanly.
