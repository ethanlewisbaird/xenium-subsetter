# xenium-subsetter

Subset a 10x Genomics Xenium dataset to a user-defined region of interest (ROI), producing a self-contained directory ready for downstream segmentation tools (BIDCell, Segger, Baysor, etc.) and for creating a xeniumranger-compatible bundle for import into Xenium Explorer.

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
            │  (using segmentation outputs from your method)
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
- `morphology_focus/ch0000_dapi.ome.tif` — cropped DAPI image
- `nuclei.tif` / `dapi_resized.tif` — downsampled for BIDCell
- `affine.csv` — affine transform for BIDCell
- `subset_offsets.json` — records the `(x_offset_um, y_offset_um)` from original space

### Step 3 — Build a xeniumranger bundle

```bash
python -m xenium_subsetter.build_bundle \
    --subset-dir  /path/to/outs_subset \
    --xenium-dir  /path/to/outs \
    --output-dir  /path/to/outs_subset_bundle
```

The bundle symlinks the large `transcripts.zarr.zip` (3–4 GB) from the original `outs/` and uses the cropped morphology. Pass `--copy` to create a fully self-contained copy instead.

### Step 4 — Run your segmentation method

Use `outs_subset/` as input to segmentation tools (coordinates are in subset space).

### Step 5 — Import into Xenium Explorer

```bash
xeniumranger import-segmentation \
    --id            experiment_<method> \
    --xenium-bundle outs_subset_bundle \
    --transcript-assignment segmentation_baysor.csv \
    --viz-polygons  segmentation_polygons.geojson \
    --units         microns \
    --localcores    8 \
    --localmem      48
```

See the [xenium-benchmark](../xenium-benchmark) project for converters that produce the correct input formats.

## Coordinate system

The subset uses **subset-relative coordinates**: the top-left corner of the ROI (minus the margin) is `(0, 0)`. Downstream segmentation tools work in this space.

When importing into Xenium Explorer via xeniumranger, cell boundaries and transcript assignments must be converted back to **original Xenium coordinates** by adding `x_offset_um` and `y_offset_um` from `subset_offsets.json`. The converter scripts in [xenium-benchmark](../xenium-benchmark) handle this automatically.

**Note on morphology alignment**: The `transcripts.zarr.zip` in the bundle is from the original full-slide Xenium output and uses original coordinates. The cropped morphology image starts at pixel `(0,0)`, which Xenium Explorer maps to µm `(0,0)`. This means the morphology image will appear offset from the transcripts/cells by `(x_offset_um, y_offset_um)` in the final Explorer view. For benchmarking purposes this is generally acceptable. Full coordinate-aligned subset bundles (requiring subsetting of the transcripts zarr) are a planned future feature.

## Output directory structure

```
outs_subset/
├── transcripts.parquet          subset transcripts, subset coords
├── transcripts.csv.gz           same in CSV
├── cell_boundaries.parquet      10x original cell boundaries for subset
├── nucleus_boundaries.parquet
├── morphology_focus/
│   └── ch0000_dapi.ome.tif      cropped DAPI (2125 Å/px)
├── nuclei.tif                   downsampled for BIDCell (5000 Å/px)
├── dapi_resized.tif             same
├── affine.csv                   BIDCell affine transform
└── subset_offsets.json          {x_offset_um, y_offset_um, ...}
```
