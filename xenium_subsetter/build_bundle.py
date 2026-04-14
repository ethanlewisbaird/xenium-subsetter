"""
Build a valid Xenium bundle for use with xeniumranger import-segmentation.

Takes an outs_subset/ directory (from subset.py) and creates a
bundle that xeniumranger accepts as --xenium-bundle input.

Key constraint: xeniumranger requires transcripts.zarr.zip in the bundle.
This file is large (~3.3 GB). Two strategies are supported:

  Strategy A (default): symlink the full-slide transcripts.zarr.zip from the
    original outs/ directory. Saves disk space. Requires the original outs/ to
    remain accessible while running xeniumranger.

  Strategy B (--copy): copy the zarr. Larger disk footprint but fully
    self-contained.

NOTE on coordinate alignment
----------------------------
The transcripts.zarr.zip is in original Xenium coordinate space.
The subset morphology image starts at pixel (0,0) = µm (0,0) in Xenium Explorer,
which means the morphology will appear offset from the transcripts by
(x_offset_um, y_offset_um) when viewed in Xenium Explorer.

For correct alignment, segmentation outputs (cell boundaries, viz polygons)
must be converted to ORIGINAL Xenium coordinates before running xeniumranger.
See converters in the xenium-benchmark project for this transform.

The subset_offsets.json written by subset.py records the required offsets.

Usage
-----
    python -m xenium_subsetter.build_bundle \\
        --subset-dir /path/to/outs_subset \\
        --xenium-dir /path/to/outs \\          # for the transcripts.zarr.zip
        --output-dir /path/to/outs_subset_bundle \\
        [--copy]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path


def build_xenium_bundle(
    subset_dir: str | Path,
    xenium_dir: str | Path,
    output_dir: str | Path,
    copy_zarr: bool = False,
) -> Path:
    """
    Build a xeniumranger-compatible bundle from a subset directory.

    Parameters
    ----------
    subset_dir  : Path to outs_subset/ (created by subset.py)
    xenium_dir  : Path to original Xenium outs/ (for transcripts.zarr.zip,
                  cells.zarr.zip, analysis_summary.html, gene_panel.json)
    output_dir  : Where to write the bundle
    copy_zarr   : If True, copy zarr files (self-contained but uses more disk)

    Returns
    -------
    Path to the created bundle directory
    """
    subset_dir = Path(subset_dir)
    xenium_dir = Path(xenium_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def link_or_copy(src: Path, dst: Path, always_copy: bool = False) -> None:
        if dst.exists():
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        if always_copy or copy_zarr:
            shutil.copy2(src, dst)
        else:
            try:
                os.link(src, dst)  # hardlink (zero extra disk)
            except OSError:
                os.symlink(src.resolve(), dst)  # symlink as fallback

    def link_dir(src: Path, dst: Path, always_copy: bool = False) -> None:
        if dst.exists():
            return
        if always_copy or copy_zarr:
            shutil.copytree(src, dst)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.symlink(src.resolve(), dst)
            except OSError:
                shutil.copytree(src, dst)

    # ── Large files from original outs/ (link by default) ─────────────────────
    large_files = [
        "transcripts.zarr.zip",   # required by xeniumranger
        "cells.zarr.zip",         # required (original 10x segmentation)
        "analysis_summary.html",  # required for payload extraction
    ]
    for fname in large_files:
        src = xenium_dir / fname
        if src.exists():
            link_or_copy(src, output_dir / fname)
            print(f"  Linked {fname} ({src.stat().st_size / 1e9:.2f} GB)")
        else:
            print(f"  WARNING: {fname} not found in {xenium_dir}")

    # Large dirs
    for dname in ["aux_outputs"]:
        src = xenium_dir / dname
        if src.exists():
            link_dir(src, output_dir / dname)
            print(f"  Linked {dname}/")

    # ── Small shared files (always copy) ─────────────────────────────────────
    small_files = ["gene_panel.json", "nucleus_boundaries.parquet",
                   "nucleus_boundaries.csv.gz"]
    for fname in small_files:
        src = xenium_dir / fname
        if src.exists():
            link_or_copy(src, output_dir / fname, always_copy=True)

    # ── Subset-specific files (copy from subset_dir) ─────────────────────────
    # Morphology: use the cropped subset image (NB: coordinate offset applies,
    # see module docstring)
    morph_src = subset_dir / "morphology_focus" / "ch0000_dapi.ome.tif"
    morph_dst = output_dir / "morphology_focus" / "ch0000_dapi.ome.tif"
    if morph_src.exists():
        link_or_copy(morph_src, morph_dst, always_copy=False)
        print(f"  Linked subset morphology ({morph_src.stat().st_size / 1e9:.2f} GB)")
    else:
        # Fall back to full-slide morphology
        full_morph = xenium_dir / "morphology_focus" / "ch0000_dapi.ome.tif"
        if full_morph.exists():
            link_or_copy(full_morph, morph_dst)
            print(f"  WARNING: subset morphology not found; using full-slide morphology")

    # Nucleus boundaries from subset
    for fname in ["nucleus_boundaries.parquet"]:
        src = subset_dir / fname
        if src.exists():
            link_or_copy(src, output_dir / fname, always_copy=True)

    # ── Build experiment.xenium ───────────────────────────────────────────────
    src_xenium = xenium_dir / "experiment.xenium"
    if not src_xenium.exists():
        raise FileNotFoundError(f"experiment.xenium not found in {xenium_dir}")

    import json as _json
    xenium_data = _json.loads(src_xenium.read_text())

    # Update morphology path to use subset image
    xenium_data["images"]["morphology_focus_filepath"] = "morphology_focus/ch0000_dapi.ome.tif"

    # Ensure paths are relative and files exist
    xenium_data["xenium_explorer_files"]["transcripts_zarr_filepath"] = "transcripts.zarr.zip"

    # Write updated experiment.xenium
    out_xenium = output_dir / "experiment.xenium"
    out_xenium.write_text(_json.dumps(xenium_data, indent=2))
    print(f"  Written experiment.xenium")

    # Copy subset offsets for reference
    offsets_src = subset_dir / "subset_offsets.json"
    if offsets_src.exists():
        shutil.copy2(offsets_src, output_dir / "subset_offsets.json")
        offsets = _json.loads(offsets_src.read_text())
        print(
            f"\n  IMPORTANT: subset offsets are "
            f"x={offsets['x_offset_um']:.2f} µm, y={offsets['y_offset_um']:.2f} µm\n"
            f"  Segmentation outputs must be converted to original Xenium coordinates\n"
            f"  before running xeniumranger. Use the converters in xenium-benchmark."
        )

    print(f"\nBundle created at: {output_dir}")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    parser.add_argument("--subset-dir", required=True, help="outs_subset/ from subset.py")
    parser.add_argument("--xenium-dir", required=True, help="Original Xenium outs/ directory")
    parser.add_argument("--output-dir", required=True, help="Output bundle directory")
    parser.add_argument("--copy", action="store_true",
                        help="Copy zarr files instead of symlinking (larger disk usage)")
    args = parser.parse_args()

    build_xenium_bundle(
        subset_dir=args.subset_dir,
        xenium_dir=args.xenium_dir,
        output_dir=args.output_dir,
        copy_zarr=args.copy,
    )


if __name__ == "__main__":
    main()
