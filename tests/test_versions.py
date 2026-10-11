"""The version has to agree in every file that declares it."""

from __future__ import annotations

from pathlib import Path

import tomllib

import molecular_loader

REPOSITORY = Path(__file__).parent.parent


def test_version_is_the_same_everywhere():
    with open(REPOSITORY / "pyproject.toml", "rb") as handle:
        project = tomllib.load(handle)["project"]["version"]
    with open(REPOSITORY / "blender_manifest.toml", "rb") as handle:
        manifest = tomllib.load(handle)["version"]

    assert molecular_loader.__version__ == project
    assert manifest == project
    assert ".".join(str(part) for part in molecular_loader.bl_info["version"]) == project


def test_addon_metadata_is_complete():
    for key in ("name", "author", "version", "blender", "location", "description", "category"):
        assert molecular_loader.bl_info[key], f"bl_info is missing {key}"
