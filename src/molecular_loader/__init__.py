"""Load molecules and molecular dynamics trajectories into Blender.

The package is both the library and the Blender add-on: activating it adds
``File > Import > Molecular Structure`` and a *Molecular Loader* panel to the
3D viewport sidebar, so no file has to be edited before loading a structure.

From a script or from Blender's Text Editor the same thing is one call:

>>> import molecular_loader
>>> molecular_loader.load("molecule.xyz")
>>> molecular_loader.load_from_environment()   # uses the MOLECULAR_LOADER_* variables

``ase`` (and the ``scipy`` it imports) have to be importable from Blender's
Python; see :mod:`molecular_loader.deps` and the README.
"""

from __future__ import annotations

import importlib
from typing import Any

from molecular_loader.deps import add_library_directories

bl_info = {
    "name": "Molecular Loader",
    "author": "ortin-fj",
    "version": (0, 2, 0),
    "blender": (3, 6, 0),
    "location": "File > Import, and 3D Viewport > Sidebar > Molecular Loader",
    "description": "Load molecules and trajectories as spheres, bonds and keyframes",
    "category": "Import-Export",
    "doc_url": "https://github.com/ortin-fj/molecular_loader",
}

__version__ = "0.2.0"

# The dependencies may live in the ``libs`` directory of the add-on, so make them
# importable before any submodule needs them.  This does not import Blender or ASE.
add_library_directories()

#: Public names, imported on demand so that ``import molecular_loader`` works
#: without Blender (and without ASE) present.
_LAZY_EXPORTS = {
    "ATOM_PROPERTIES": ("molecular_loader.elements", "ATOM_PROPERTIES"),
    "LoadReport": ("molecular_loader.loader", "LoadReport"),
    "LoaderSettings": ("molecular_loader.settings", "LoaderSettings"),
    "bond_geometry": ("molecular_loader.geometry", "bond_geometry"),
    "find_bonds": ("molecular_loader.bonds", "find_bonds"),
    "load": ("molecular_loader.loader", "load"),
    "load_from_environment": ("molecular_loader.loader", "load_from_environment"),
    "load_structure": ("molecular_loader.loader", "load_structure"),
    "parse_bond_pairs": ("molecular_loader.settings", "parse_bond_pairs"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    try:
        module_name, attribute = _LAZY_EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None

    value = getattr(importlib.import_module(module_name), attribute)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))


def register() -> None:
    """Entry point of the Blender add-on."""
    from molecular_loader import ui

    ui.register()


def unregister() -> None:
    """Called by Blender when the add-on is disabled."""
    from molecular_loader import ui

    ui.unregister()


if __name__ == "__main__":
    register()
