# xenium-subsetter

Subset a 10x Genomics Xenium dataset to a user-defined region of interest (ROI), producing a self-contained directory ready for downstream segmentation tools (BIDCell, Segger, Baysor, etc.) and a xeniumranger-compatible bundle for import into Xenium Explorer.

## Motivation

Full Xenium slides can be 10+ GB and cover large areas. When benchmarking segmentation methods on a specific biological region, it is much more efficient to work on a cropped subset:

- Faster tool runtimes
- Smaller disk/memory footprint
- Easier sharing and reproducibility

## Workflow

```
Xenium Explorer
    └── Export selection polygon CSVs
            │
            ▼
   xenium_subsetter.subset
            │  crops transcripts, boundaries, DAPI
            ▼
      outs_subset/         ← subset-space coordinates (origin = ROI top-left)
            │
   xenium_subsetter.build_bundle
            │  links/copies files into a valid Xenium bundle
            ▼
  outs_subset_bundle/      ← valid --xenium-bundle for xeniumranger
            │
  xeniumranger import-segmentation
            │  (via xenium-segmentation-benchmark)
            ▼
  experiment_<method>/     ← Xenium Explorer .xenium file
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Step 1 — Define the region in Xenium Explorer

1. Open your dataset in [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer)
2. Draw a selection polygon around your region of interest
3. Export the selection: **File → Export Selection → Save as CSV**
4. Repeat for multiple regions if desired; place all CSVs in one folder

### Step 2 — Create the subset

```bash
python -m xenium_subsetter.subset \
    --xenium-dir  /path/to/outs \
    --coords-dir  /path/to/selection_coordinates \
    --output-dir  /path/to/outs_subset \
    --margin      200
```

This produces `outs_subset/` containing:
- `transcripts.parquet` — transcripts in subset-relative coordinates
- `transcripts.csv.gz` — same data in CSV format (for tools that need it)
- `cell_boundaries.parquet` / `nucleus_boundaries.parquet` — boundary polygons
- `morphology_focus/ch0000_dapi.ome.tif` — cropped DAPI image (with OME-XML UUID for xeniumranger 4.x compatibility)
- `nuclei.tif` / `dapi_resized.tif` — downsampled for BIDCell
- `affine.csv` — affine transform for BIDCell
- `subset_offsets.json` — records `{x_offset_um, y_offset_um}` for converting back to original coordinates

### Step 3 — Build a xeniumranger bundle

```bash
python -m xenium_subsetter.build_bundle \
    --subset-dir  /path/to/outs_subset \
    --xenium-dir  /path/to/outs \
    --output-dir  /path/to/outs_subset_bundle
```

The bundle symlinks the large `transcripts.zarr.zip` (3–4 GB) from the original `outs/` and uses the cropped morphology image. Also links `morphology.ome.tif` and `metrics_summary.csv` required by xeniumranger ≥ 4.0. Pass `--copy` for a fully self-contained copy.

### Step 4 — Run your segmentation method

Use `outs_subset/` as input to segmentation tools. Coordinates are in subset-relative space (origin = ROI top-left minus margin).

See [xenium-segmentation-benchmark](https://github.com/ethanlewisbaird/xenium-segmentation-benchmark) for a complete pipeline wrapping Segger and BIDCell.

### Step 5 — Import into Xenium Explorer

Segmentation outputs must be converted back to original Xenium coordinates before import. The `to_xenium.py` scripts in xenium-segmentation-benchmark handle this using `subset_offsets.json`.

```bash
xeniumranger import-segmentation \
    --id            experiment_<method> \
    --xenium-bundle /path/to/outs_subset_bundle \
    --cells         /path/to/cells.geojson \
    --units         microns \
    --localcores    8 \
    --localmem      48
```

## Coordinate system

The subset uses **subset-relative coordinates**: the top-left corner of the ROI (minus the margin) is `(0, 0)`. All segmentation tools work in this space.

When importing into Xenium Explorer via xeniumranger, cell boundaries and transcript assignments must be converted back to **original Xenium coordinates** by adding `x_offset_um` and `y_offset_um` from `subset_offsets.json`.

**Note on morphology alignment**: `transcripts.zarr.zip` in the bundle uses original Xenium coordinates, while the cropped morphology image starts at pixel `(0, 0)` which Xenium Explorer maps to µm `(0, 0)`. The morphology will appear offset from transcripts/cells by the subset offset when viewed in Explorer. For benchmarking segmentation quality this is generally acceptable.

## Output directory structure

```
outs_subset/
├── transcripts.parquet          subset transcripts, subset coords
├── transcripts.csv.gz           same in CSV
├── cell_boundaries.parquet      original 10x cell boundaries clipped to subset
├── nucleus_boundaries.parquet
├── morphology_focus/
│   └── ch0000_dapi.ome.tif      cropped DAPI (0.2125 µm/px, OME-TIFF with UUID)
├── nuclei.tif                   downsampled for BIDCell (0.5 µm/px)
├── dapi_resized.tif             same
├── affine.csv                   BIDCell affine transform
└── subset_offsets.json          {x_offset_um, y_offset_um, ...}
```
