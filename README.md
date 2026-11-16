# lightsuite-neuroglancer

Open LightSuite atlas-registration results in [Neuroglancer](https://github.com/google/neuroglancer):
your sample channels with the warped atlas on top, with region names, colours,
the atlas hierarchy, per-region statistics and 3D region meshes.

## Usage

Needs [uv](https://docs.astral.sh/uv/) (Linux, Windows or macOS). Python and the
dependencies are downloaded automatically on first run (internet needed once).
Tested on Linux and Windows 11. Also needed: a browser with WebGL 2 (recent
Chrome, Firefox or Edge), 3–4 GB free disk per brain, and probably 8 GB or more
of RAM for the first conversion of a 20 µm brain.

```bash
cd lightsuite-neuroglancer
uv run view_in_neuroglancer.py                       # opens a folder picker
uv run view_in_neuroglancer.py /path/to/results       # or give the folder
```

Or install it as a command (works on Windows, Linux and macOS):

```bash
uv tool install --python 3.12 git+https://github.com/YCabara/lightsuite_neuroglancer.git
lightsuite-neuroglancer D:\brain8_analysed\8_DNF_atlased
```

`pip install .` into a Python 3.10–3.13 environment also works. Newer Pythons lack
prebuilt packages for some dependencies and would need a C++ compiler.

**Windows, no setup:** download the repository (ZIP), extract it and double-click
`run_viewer.bat`, or drop a results folder on it. If uv is missing, it installs it,
then it uses Python 3.12 and the tested package versions automatically.

Documentation:

* [docs/installation.html](docs/installation.html): installation, options, troubleshooting
* [docs/user-guide.html](docs/user-guide.html): how to check a registration, find regions, read statistics
* [docs/troubleshooting.html](docs/troubleshooting.html): checklist for a proper run, region names, error messages
* [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): code structure, dependencies, maintenance (for developers)

The folder can be the registration output itself (the one containing
`volume_registered/sample_space/sample_space_manifest.json`) or any folder
above it. If several results are found you are asked to pick one.

The first run converts the data (a few minutes for a whole brain, mostly
building meshes) into `<results>/neuroglancer/` (about 3–4 GB for a brain at
20 µm). Later runs reuse it and start immediately. It is rebuilt automatically
when the registration files or this tool change. A local viewer then opens in
your browser. Leave the terminal running while you use it, and press Ctrl+C to stop.

Options:

| option | |
|---|---|
| `--atlas-dir DIR` | folder with the atlas ontology CSV if it is not found automatically |
| `--out DIR` | write the converted data somewhere else |
| `--no-meshes` | skip 3D meshes (much faster conversion) |
| `--zoom DIR` | also show a LightSuite multires ROI (e.g. a 3.2x stack) at native resolution, see [Zoom stacks](#zoom-stacks); repeatable |
| `--reexport` | force re-conversion (atlas and zoom data) |
| `--layout 3d` | start layout: `4panel` (default), `xy`, `yz`, `xz`, `3d` |
| `--port N` | port of the local data server (default: any free port) |
| `--no-browser` | only print the URL |

## What you see

| layer | |
|---|---|
| `ch1`, `ch2`, … | sample channels on the registration grid (ch2+ hidden, coloured) |
| `template` | atlas template warped onto the sample (hidden) |
| `boundaries` | atlas region outlines |
| one layer per hierarchy level | e.g. Allen `substructure`, `structure`, `division`, `category` |
| extra label volumes | e.g. spinal cord `segments`, `division_labels` (hidden) |
| `region meshes` | 3D surface of every region at every level (whole brain/cord shown) |

* Hover to see the region name at every level. The **Selection** panel
  (bottom right) lists all statistics of the region under the mouse.
* **Double-click** a region: only that region stays shown, its 3D mesh is added
  to `region meshes`, and its statistics are pinned in the Selection panel, for
  the region and each region containing it (e.g. CP → STR → CNU). Double-click it
  again to show all regions. Double-click the background to unpin.
  (Ctrl+right-click pins without isolating the region.)
* In a layer's **Seg.** tab, search by acronym, name, ID, or `#tag`
  (tags are the division/category a region belongs to, e.g. `#Isocortex`).
  Click a region to highlight it. In `region meshes` this shows its 3D surface.
* Per-region statistics also show under "numerical properties" in the Seg.
  tab, where you can add them as columns, and sort and filter by them.

### Side panel: region search and channel LUTs

The browser opens on `http://127.0.0.1:<port>/lightsuite/`: Neuroglancer with a
panel on the left (hide it with ‹). The plain Neuroglancer URL is printed too.

* **Find a region**: type a name, acronym or ID, e.g. `locus coeruleus`, `CA1`,
  `147`. Matching is forgiving: case, accents, British spellings
  (coeruleus = ceruleus) and small typos are tolerated; weaker matches are marked
  "similar". **Enter** (or a click) shows the first / chosen match: that region
  alone is highlighted on its hierarchy level, the view jumps to it and zooms to
  fit it, its 3D mesh is shown and its statistics (and those of the regions
  containing it) are pinned in the Selection panel. **Shift+click** /
  **Shift+Enter** adds a region to the highlight. **Show all regions** restores
  the previous view.
* **Channels**: per channel (and per zoom-stack channel) a visibility checkbox,
  a LUT — Grays, Green, Red, Blue, Magenta, Cyan, Yellow, Fire, Inferno, Magma,
  Plasma, Viridis, Turbo, Jet, HiLo (blue = at/below the range minimum,
  red = at/above the maximum, to check saturation) — and **invert**. Changing the
  LUT keeps the contrast range set in the layer's Rendering tab. Channels 2+ are
  blended additively, so several channels can be shown together.

## Where region names come from

Tried in this order:

1. **Spinal cord (Fiederling) atlas**: `Atlas_Regions.csv` next to the
   manifest's `segments_csv`.
2. **Allen CCF**: `parcellation_to_parcellation_term_membership.csv` in
   `--atlas-dir`, the config's `atlas.atlas_dir` (or their subfolders) or `~/atlases/*`.
3. **Perens / Gubra LSFM**: `ARA2_annotation_info.csv` in the same folders.
4. **BrainGlobe atlas**: `structures.json` of a downloaded atlas in
   `~/.brainglobe`, `/home/*/.brainglobe` or `$BRAINGLOBE_DIR`,
   matched by the config's `brainglobe_name` (e.g. Perens LSFM, Princeton).
5. **The sample's own** `stats/region_stats.csv` (region / structure / division).
6. Bare label numbers.

## Region statistics

All of these are added when present (hover over a statistic's name in the
Selection panel to see where it comes from):

| source | names |
|---|---|
| `registration_analysis/region_metrics.csv` | every numeric column, e.g. `voxel_count`, `volume_mm3`, `mean_intensity`, `positive_percent` (every hierarchy level) |
| `stats/region_stats_sample.csv` (sample space) | `median_intensity_left`, `std_right`, `volume_mm3_left`, … (finest level) |
| `stats/region_stats.csv` (atlas space) | the same, prefixed `atlas_`, e.g. `atlas_volume_mm3_left` = volume of the atlas region |

Several channels add `_ch1`, `_ch2`, … Columns with one value per channel
(threshold, background) and empty columns are left out. Regions without a value
show nothing (they are not shown as 0). The spinal cord per-segment tables are
not shown yet.

Labels that are not in the atlas are listed as `label N (not in atlas)`.

## Orientation

Volumes are shown on the pipeline's registration grid, the same way
LightSuite's own QC overlays show them. Channels get
`permute_sample_to_atlas` applied. Atlas volumes get it too, unless the manifest
says `atlas_volumes_permuted: true`.

## Zoom stacks

```bash
uv run view_in_neuroglancer.py /path/to/atlased --zoom /path/to/multires_results
```

`--zoom` takes a LightSuite multires results folder (the one with
`multires_regopts.json`). Each ROI channel of its pair manifest gets a layer,
e.g. `greencells 488` and `greencells 561`, placed over the whole brain and below
the atlas outlines and labels. The view starts centred on the first zoom. With the default 4panel layout the 3D panel is replaced by an xy view zoomed on
the (first) ROI: the zoom at native resolution with the atlas outlines and labels,
no region meshes. It follows the shared position but keeps its own zoom level.

* **Nothing is registered or resampled again.** The raw ROI TIFFs (the
  `roi.volume_path` of each channel) are converted as they are into a
  multiscale precomputed volume. `s0` is the native voxel size (e.g.
  2.03 x 2.03 x 3 µm), and each further level halves it. LightSuite's resampled
  `*_registered_to_*.tif` files are not used.
* **They are placed with one affine** built from parameters the multires run
  already saved: the voxel size, stage position and orientation in the pair manifest
  (`converted/*_pair.json`), the landmark fit `roi_to_overview_tform` in
  `multires_regopts.json`, and the elastix `TransformParameters.*.txt` chain
  (Translation / Euler / Affine). These are followed by the atlas run's
  registration-grid downsampling (`regopts.json`) and `permute_sample_to_atlas`.
  Neuroglancer applies the matrix when it draws the layer.
  A B-spline elastix stage is refused, because it cannot be shown as one matrix.
* The multires overview stack must be the stack the atlas registration used
  (same size and voxel size), otherwise the tool stops.
* Converted data goes to `<multires_results>/neuroglancer/` (about 1.15x the
  raw TIFF size, uncompressed) and is reused until the TIFFs change. The
  matrix is recomputed at every start and also written to `zoom_meta.json`.
* The atlas labels are still on their 20 µm grid, so region edges look blocky when
  you are zoomed in on the stack.
