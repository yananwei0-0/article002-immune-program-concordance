#!/usr/bin/env python3
"""Execute the read-only TCGA source script with v2-safe path rendering.

The historical implementation assumes every report asset is below the data
lake root. v2 outputs must instead stay in the authorized analysis directory;
only the display helper is replaced. Statistical functions and arguments are
unchanged.
"""

from __future__ import annotations
from sci27_portability import expand_legacy_paths

import importlib.util
import sys
from pathlib import Path


SOURCE = Path(
    expand_legacy_paths("${SCI27_RESEARCH_ROOT}/outputs/manuscript17_open_science_packages_20260803/"
    "BIO-05_pan_cancer_proteogenomic_immune_evasion/repository/code/project/scripts/"
    "run_pan_cancer_proteomic_immune_evasion_v0_5_tcga_sensitivity.py")
)


def main() -> None:
    spec = importlib.util.spec_from_file_location("bio05_tcga_readonly", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(SOURCE.parent))
    spec.loader.exec_module(module)

    def safe_rel(path: Path, root: Path) -> str:
        try:
            return str(path.relative_to(root))
        except ValueError:
            return str(path)

    module.rel = safe_rel
    module.main()


if __name__ == "__main__":
    main()
