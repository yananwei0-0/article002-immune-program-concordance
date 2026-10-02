"""Portable path expansion for the SCI27 repository package."""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def _legacy(*parts: str) -> str:
    return str(Path("/").joinpath(*parts))


def _configured(name: str, fallback: Path) -> str:
    return str(Path(os.environ.get(name, str(fallback))).expanduser().resolve())


def expand_legacy_paths(value: str) -> str:
    """Map legacy absolute prefixes to configured release roots."""
    external = REPO_ROOT / "external_data"
    token_mapping = {
        "${SCI27_REPO_ROOT}": _configured("SCI27_REPO_ROOT", REPO_ROOT),
        "${SCI27_RESEARCH_ROOT}": _configured("SCI27_RESEARCH_ROOT", REPO_ROOT),
        "${SCI27_DATA_ROOT}": _configured("SCI27_DATA_ROOT", external),
        "${SCI27_BIOBANK_ROOT}": _configured("SCI27_BIOBANK_ROOT", external / "biobank"),
        "${SCI27_TMP_ROOT}": _configured("SCI27_TMP_ROOT", Path(tempfile.gettempdir())),
    }
    mapping = {
        _legacy("Users", "yanan", "Desktop", "research"): _configured("SCI27_RESEARCH_ROOT", REPO_ROOT),
        _legacy("Users", "yanan", "Desktop", "bioinformatics_analysis"): _configured("SCI27_DATA_ROOT", external),
        _legacy("Volumes", "Biobank"): _configured("SCI27_BIOBANK_ROOT", external / "biobank"),
        _legacy("private", "tmp"): _configured("SCI27_TMP_ROOT", Path(tempfile.gettempdir())),
        _legacy("tmp"): _configured("SCI27_TMP_ROOT", Path(tempfile.gettempdir())),
    }
    result = str(value)
    for marker, replacement in token_mapping.items():
        result = result.replace(marker, replacement)
    pattern = re.compile("|".join(re.escape(item) for item in sorted(mapping, key=len, reverse=True)))
    result = pattern.sub(lambda match: mapping[match.group(0)], result)
    return os.path.expandvars(result)


def expand_environment_tokens(value):
    """Recursively expand ${SCI27_*} tokens in rendered configuration data."""
    if isinstance(value, str):
        return os.path.expandvars(value)
    if isinstance(value, list):
        return [expand_environment_tokens(item) for item in value]
    if isinstance(value, dict):
        return {key: expand_environment_tokens(item) for key, item in value.items()}
    return value
