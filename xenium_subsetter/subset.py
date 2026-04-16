"""
Subset a 10x Xenium dataset to a user-defined region of interest.

Reads selection polygons exported from Xenium Explorer, computes a bounding
box with a configurable margin, and produces a self-contained output directory
with coordinate-shifted files ready for BIDCell, Segger, and other tools.

The output uses subset-relative coordinates (origin at the top-left corner of
the bounding box), so downstream tools see a "normal" dataset starting at (0,0).

Usage
-----
    python -m xenium_subsetter.subset \\
        --xenium-dir /path/to/outs \\
        --coords-dir /path/to/selection_coordinates_xenium_explorer \\
        --output-dir /path/to/outs_subset \\
        [--margin 200]

Outputs
-------
    outs_subset/
    ├── transcripts.parquet       (subset, coords shifted to 0-origin)
    ├── transcripts.csv.gz        (same, CSV format for tools that need it)
    ├── cell_boundaries.parquet   (original 10x cell boundaries for subset)
    ├── nucleus_boundaries.parquet
    ├── morphology_focus/
    │   └── ch0000_dapi.ome.tif   (cropped DAPI image for subset region)
    ├── nuclei.tif                (resized DAPI for BIDCell nucleus segmentation)
    ├── dapi_resized.tif          (same as nuclei.tif, alternative name)
    └── affine.csv                (affine transform for BIDCell)
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import tifffile


def read_selection_polygons(coords_dir: str | Path) -> tuple[list[float], list[float]]:
    """Read x,y coordinates from all Xenium Explorer selection CSVs."""
    xs, ys = [], []
    for fp in sorted(glob.glob(os.path.join(str(coords_dir), "*.csv"))):
        for line in open(fp, encoding="latin-1").readlines()[3:]:
            line = line.strip().strip('"')
            if "," in line:
                try:
                    x, y = line.split(",")
                    xs.append(float(x))
                    ys.append(float(y))
                except ValueError:
                    continue
    if not xs:
        raise ValueError(f"No coordinate data found in {coords_dir}")
    return xs, ys


def compute_bbox(
    xs: list[float],
    ys: list[float],
    margin_um: float = 200.0,
    pixel_size_um: float = 0.2125,
) -> dict:
    """Compute bounding box from selection polygon coordinates."""
    x_min_um = max(0.0, min(xs) - margin_um)
    x_max_um = max(xs) + margin_um
    y_min_um = max(0.0, min(ys) - margin_um)
    y_max_um = max(ys) + margin_um

    return {
        "x_min_um": x_min_um,
        "x_max_um": x_max_um,
        "y_min_um": y_min_um,
        "y_max_um": y_max_um,
        "row_start": int(y_min_um / pixel_size_um),
        "row_end": int(y_max_um / pixel_size_um) + 1,
        "col_start": int(x_min_um / pixel_size_um),
        "col_end": int(x_max_um / pixel_size_um) + 1,
        "pixel_size_um": pixel_size_um,
    }


def subset_transcripts(
    input_dir: Path,
    output_dir: Path,
    bbox: dict,
) -> None:
    """Subset transcripts.parquet (and .csv.gz) to the bounding box."""
    x_min, x_max = bbox["x_min_um"], bbox["x_max_um"]
    y_min, y_max = bbox["y_min_um"], bbox["y_max_um"]

    # Load and filter
    tx_path = input_dir / "transcripts.parquet"
    if not tx_path.exists():
        raise FileNotFoundError(f"transcripts.parquet not found in {input_dir}")

    print(f"  Loading transcripts from {tx_path}...")
    df = pq.read_table(tx_path).to_pandas()
    mask = (
        (df["x_location"] >= x_min) & (df["x_location"] <= x_max) &
        (df["y_location"] >= y_min) & (df["y_location"] <= y_max)
    )
    sub = df[mask].copy()

    # Shift coordinates to subset origin
    sub["x_location"] -= x_min
    sub["y_location"] -= y_min
    if "x_local_px" in sub.columns:
        sub["x_local_px"] -= bbox["col_start"]
        sub["y_local_px"] -= bbox["row_start"]

    print(f"  Retained {len(sub):,} / {len(df):,} transcripts")

    out_pq = output_dir / "transcripts.parquet"
    pq.write_table(pa.Table.from_pandas(sub, preserve_index=False), out_pq)
    print(f"  Written to {out_pq}")

    # Also write CSV.gz for tools that need it
    out_csv = output_dir / "transcripts.csv.gz"
    sub.to_csv(out_csv, index=False, compression="gzip")
    print(f"  Written to {out_csv}")


def subset_boundaries(
    input_dir: Path,
    output_dir: Path,
    bbox: dict,
    filename: str,
) -> None:
    """Subset a boundary parquet file to the bounding box."""
    src = input_dir / filename
    if not src.exists():
        print(f"  WARNING: {filename} not found, skipping")
        return

    x_min, x_max = bbox["x_min_um"], bbox["x_max_um"]
    y_min, y_max = bbox["y_min_um"], bbox["y_max_um"]

    df = pq.read_table(src).to_pandas()
    mask = (
        (df["vertex_x"] >= x_min) & (df["vertex_x"] <= x_max) &
        (df["vertex_y"] >= y_min) & (df["vertex_y"] <= y_max)
    )
    sub = df[mask].copy()
    sub["vertex_x"] -= x_min
    sub["vertex_y"] -= y_min

    dst = output_dir / filename
    pq.write_table(pa.Table.from_pandas(sub, preserve_index=False), dst)
    print(f"  Written {filename}: {len(sub):,} rows")


def crop_morphology(
    input_dir: Path,
    output_dir: Path,
    bbox: dict,
) -> None:
    """Crop the DAPI morphology focus image to the bounding box."""
    src = input_dir / "morphology_focus" / "ch0000_dapi.ome.tif"
    if not src.exists():
        print(f"  WARNING: morphology_focus/ch0000_dapi.ome.tif not found, skipping")
        return

    r0, r1 = bbox["row_start"], bbox["row_end"]
    c0, c1 = bbox["col_start"], bbox["col_end"]

    print(f"  Cropping morphology: rows {r0}-{r1}, cols {c0}-{c1}...")
    img = tifffile.imread(src)
    # Handle multi-dimensional images
    if img.ndim == 2:
        crop = img[r0:r1, c0:c1]
    elif img.ndim == 3:
        crop = img[:, r0:r1, c0:c1]
    else:
        crop = img[..., r0:r1, c0:c1]

    out_dir = output_dir / "morphology_focus"
    out_dir.mkdir(exist_ok=True)
    dst = out_dir / "ch0000_dapi.ome.tif"

    # xeniumranger 4.x requires a <UUID> element inside <TiffData> in the
    # OME-XML, which plain tifffile.imwrite does not add.  Build a minimal
    # OME-XML description with a self-referencing UUID so the bundle passes
    # xeniumranger's morphology validation.
    import uuid as _uuid
    ome_uuid = f"urn:uuid:{_uuid.uuid4()}"
    h, w = (crop.shape[-2], crop.shape[-1]) if crop.ndim >= 2 else crop.shape
    ome_xml = (
        f'<?xml version="1.0" encoding="UTF-8"?>'
        f'<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06" '
        f'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
        f'xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 '
        f'http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd" '
        f'UUID="{ome_uuid}">'
        f'<Image ID="Image:0" Name="Image0"><Pixels ID="Pixels:0" '
        f'DimensionOrder="XYCZT" Type="uint16" '
        f'SizeX="{w}" SizeY="{h}" SizeC="1" SizeZ="1" SizeT="1" '
        f'PhysicalSizeX="0.2125" PhysicalSizeXUnit="µm" '
        f'PhysicalSizeY="0.2125" PhysicalSizeYUnit="µm">'
        f'<Channel ID="Channel:0:0" SamplesPerPixel="1"><LightPath/></Channel>'
        f'<TiffData IFD="0" PlaneCount="1">'
        f'<UUID FileName="ch0000_dapi.ome.tif">{ome_uuid}</UUID>'
        f'</TiffData></Pixels></Image></OME>'
    )
    tifffile.imwrite(str(dst), crop, description=ome_xml)
    print(f"  Written to {dst} ({crop.shape})")

    # Also write a resized version for BIDCell
    _write_bidcell_dapi(crop, output_dir, bbox)


def _write_bidcell_dapi(
    crop: np.ndarray,
    output_dir: Path,
    bbox: dict,
    target_pixel_size_um: float = 0.5,
) -> None:
    """Write a downsampled DAPI image for BIDCell (nuclei.tif, dapi_resized.tif)."""
    from skimage.transform import resize as sk_resize

    scale = bbox["pixel_size_um"] / target_pixel_size_um
    if crop.ndim == 2:
        h, w = crop.shape
        new_h, new_w = int(h * scale), int(w * scale)
        resized = sk_resize(crop, (new_h, new_w), preserve_range=True).astype(np.uint16)
    else:
        # 3D: take first channel / z
        plane = crop[0] if crop.ndim == 3 else crop
        h, w = plane.shape
        new_h, new_w = int(h * scale), int(w * scale)
        resized = sk_resize(plane, (new_h, new_w), preserve_range=True).astype(np.uint16)

    for fname in ["nuclei.tif", "dapi_resized.tif"]:
        tifffile.imwrite(str(output_dir / fname), resized)
    print(f"  Written BIDCell DAPI images: {resized.shape} @ {target_pixel_size_um} µm/px")


def write_affine(output_dir: Path, bbox: dict, target_pixel_size_um: float = 0.5) -> None:
    """Write affine.csv for BIDCell coordinate transform."""
    scale = target_pixel_size_um / bbox["pixel_size_um"]
    # Affine: [scale_x, 0, 0, scale_y] (no rotation/shear for axis-aligned crop)
    affine = np.array([[scale, 0, 0, scale]], dtype=np.float64)
    pd.DataFrame(affine).to_csv(output_dir / "affine.csv", header=False, index=False)
    print(f"  Written affine.csv (scale={scale:.4f})")


def subset_xenium(
    xenium_dir: str | Path,
    output_dir: str | Path,
    coords_dir: str | Path,
    margin_um: float = 200.0,
    pixel_size_um: float = 0.2125,
) -> dict:
    """
    Full subsetting pipeline.

    Parameters
    ----------
    xenium_dir   : Path to original Xenium outs/ directory
    output_dir   : Path to write the subset bundle
    coords_dir   : Path to folder containing Xenium Explorer selection CSVs
    margin_um    : Extra margin (µm) around the selection bounding box
    pixel_size_um: Pixel size of the morphology image (default 0.2125 µm/px)

    Returns
    -------
    dict with 'bbox' and 'offsets' keys describing the crop
    """
    xenium_dir = Path(xenium_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Reading selection polygons...")
    xs, ys = read_selection_polygons(coords_dir)
    bbox = compute_bbox(xs, ys, margin_um=margin_um, pixel_size_um=pixel_size_um)

    print(
        f"Bounding box (µm): "
        f"x=[{bbox['x_min_um']:.1f}, {bbox['x_max_um']:.1f}]  "
        f"y=[{bbox['y_min_um']:.1f}, {bbox['y_max_um']:.1f}]"
    )
    print(
        f"Pixel crop: rows {bbox['row_start']}-{bbox['row_end']}, "
        f"cols {bbox['col_start']}-{bbox['col_end']}"
    )

    print("\nSubsetting transcripts...")
    subset_transcripts(xenium_dir, output_dir, bbox)

    print("\nSubsetting cell boundaries...")
    subset_boundaries(xenium_dir, output_dir, bbox, "cell_boundaries.parquet")
    subset_boundaries(xenium_dir, output_dir, bbox, "nucleus_boundaries.parquet")

    print("\nCropping morphology image...")
    crop_morphology(xenium_dir, output_dir, bbox)

    print("\nWriting affine transform...")
    write_affine(output_dir, bbox)

    # Record offsets for downstream tools (so they can convert back to original coords)
    import json
    offsets = {
        "x_offset_um": bbox["x_min_um"],
        "y_offset_um": bbox["y_min_um"],
        "pixel_size_um": pixel_size_um,
        "col_offset_px": bbox["col_start"],
        "row_offset_px": bbox["row_start"],
    }
    (output_dir / "subset_offsets.json").write_text(json.dumps(offsets, indent=2))
    print(f"\nSubset offsets written to {output_dir / 'subset_offsets.json'}")
    print("Done.")
    return {"bbox": bbox, "offsets": offsets}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    parser.add_argument("--xenium-dir", required=True, help="Path to Xenium outs/ directory")
    parser.add_argument("--coords-dir", required=True,
                        help="Folder with Xenium Explorer selection CSVs")
    parser.add_argument("--output-dir", required=True, help="Output directory for subset")
    parser.add_argument("--margin", type=float, default=200.0,
                        help="Extra margin around selection (µm, default 200)")
    parser.add_argument("--pixel-size", type=float, default=0.2125,
                        help="Morphology pixel size in µm (default 0.2125)")
    args = parser.parse_args()

    subset_xenium(
        xenium_dir=args.xenium_dir,
        output_dir=args.output_dir,
        coords_dir=args.coords_dir,
        margin_um=args.margin,
        pixel_size_um=args.pixel_size,
    )


if __name__ == "__main__":
    main()
