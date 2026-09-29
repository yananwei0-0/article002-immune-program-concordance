#!/usr/bin/env python3
"""Build compact contact sheets from the regenerated 300-dpi PNG figures."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageOps


ROOT = Path(__file__).resolve().parent
PNG_DIR = ROOT / "generated" / "PNG_300dpi"
PREVIEW_DIR = ROOT / "previews"

MAIN = [
    "Fig1_CPTAC_primary_concordance",
    "Fig2_heterogeneity_and_scoring_robustness",
    "Fig3_TCGA_and_APOLLO_triangulation",
    "Fig4_single_cell_and_spatial_compartments",
]
SUPPLEMENT = [
    "FigS1_complete_pair_counts",
    "FigS2_residual_context_support",
    "FigS3_TCGA_full_correlation_atlas",
    "FigS4_TCGA_immune_subtypes",
    "FigS5_source_resolved_single_cell",
    "FigS6_threshold_and_influence_sensitivities",
]


def contact_sheet(stems: list[str], columns: int, output: Path) -> None:
    tile_width, tile_height = 780, 620
    label_height, margin = 42, 24
    rows = (len(stems) + columns - 1) // columns
    canvas = Image.new(
        "RGB",
        (columns * tile_width + (columns + 1) * margin, rows * (tile_height + label_height) + (rows + 1) * margin),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    for index, stem in enumerate(stems):
        source = PNG_DIR / f"{stem}.png"
        with Image.open(source) as raw:
            image = ImageOps.contain(raw.convert("RGB"), (tile_width, tile_height))
        row, column = divmod(index, columns)
        left = margin + column * (tile_width + margin)
        top = margin + row * (tile_height + label_height + margin)
        x = left + (tile_width - image.width) // 2
        y = top + (tile_height - image.height) // 2
        canvas.paste(image, (x, y))
        draw.text((left, top + tile_height + 8), stem, fill="#1F2937")
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, format="PNG", dpi=(150, 150), optimize=True)


def main() -> None:
    contact_sheet(MAIN, 2, PREVIEW_DIR / "main_figures_contact.png")
    contact_sheet(SUPPLEMENT, 2, PREVIEW_DIR / "supp_figures_contact.png")
    print(f"Updated contact sheets in {PREVIEW_DIR}")


if __name__ == "__main__":
    main()
