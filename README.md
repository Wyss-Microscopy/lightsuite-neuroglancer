# lightsuite-neuroglancer

Open LightSuite atlas-registration results in [Neuroglancer](https://github.com/google/neuroglancer):
your sample channels with the warped atlas on top, with region names, colours,
the atlas hierarchy, per-region statistics and 3D region meshes.

## Usage

Needs [uv](https://docs.astral.sh/uv/) (Linux, Windows or macOS). Python and the
dependencies are downloaded automatically on first run (internet needed once).
Only Linux has been tested. Also needed: a browser with WebGL 2 (recent
Chrome, Firefox or Edge), 3–4 GB free disk per brain, and probably 8 GB or more
of RAM for the first conversion of a 20 µm brain.

```bash
cd /home/yosi/lightsuite-neuroglancer
uv run view_in_neuroglancer.py                       # opens a folder picker
uv run view_in_neuroglancer.py /path/to/results       # or give the folder
```

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
| `--reexport` | force re-conversion |
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

## Where region names come from

Tried in this order:

1. **Spinal cord (Fiederling) atlas**: `Atlas_Regions.csv` next to the
   manifest's `segments_csv`.
2. **Allen CCF**: `parcellation_to_parcellation_term_membership.csv` in the
   config's `atlas.atlas_dir` or `~/atlases/*`.
3. **BrainGlobe atlas**: `structures.json` of a downloaded atlas in
   `~/.brainglobe`, `/home/*/.brainglobe` or `$BRAINGLOBE_DIR`,
   matched by the config's `brainglobe_name` (e.g. Perens LSFM, Princeton).
4. **The sample's own** `stats/region_stats.csv` (region / structure / division).
5. Bare label numbers.

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
