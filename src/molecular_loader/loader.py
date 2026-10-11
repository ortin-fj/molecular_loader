"""Build the Blender scene of a structure: the entry point of the add-on and scripts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import bpy
import numpy as np

from molecular_loader.animation import (
    extend_scene_range,
    insert_atom_keyframes,
    insert_bond_keyframes,
)
from molecular_loader.bonds import find_bonds
from molecular_loader.builder import SceneBuilder
from molecular_loader.elements import get_atom_properties
from molecular_loader.geometry import BondHalf, bond_geometry
from molecular_loader.materials import MaterialCache
from molecular_loader.settings import LoaderSettings
from molecular_loader.structure import read_snapshots, read_structure

RENDER_ENGINE = "CYCLES"


@dataclass(frozen=True)
class LoadReport:
    """What one load produced, used for the messages of the add-on and scripts."""

    filepath: Path
    collection: str
    atoms: int
    bond_halves: int
    snapshots: int
    last_frame: int

    def summary(self) -> str:
        text = (
            f"Loaded {self.atoms} atoms and {self.bond_halves} bond halves from "
            f"{self.filepath.name} into the collection '{self.collection}'"
        )
        if self.snapshots:
            return (
                f"{text}, with {self.snapshots} keyframed snapshots "
                f"(last frame {self.last_frame})."
            )
        return f"{text}; the file has a single geometry, so no keyframes were added."


def load_structure(settings: LoaderSettings, scene: bpy.types.Scene | None = None) -> LoadReport:
    """Build every object of ``settings`` and return a report of what was created."""
    scene = scene or bpy.context.scene
    scene.render.engine = RENDER_ENGINE

    structure = read_structure(settings.filepath)
    radii = np.array([get_atom_properties(symbol)["radius"] for symbol in structure.symbols])
    materials = MaterialCache()
    builder = SceneBuilder(settings.collection, structure.symbols, materials)

    atom_objects = [
        builder.add_atom(index, structure.positions[index], symbol)
        for index, symbol in enumerate(structure.symbols)
    ]

    bond_objects: list[tuple[BondHalf, bpy.types.Object]] = []
    for half in find_bonds(structure.connectivity, structure.symbols, settings.explicit_bonds):
        geometry = bond_geometry(
            structure.positions, half.owner, half.other, radii, half.same_element
        )
        if geometry is None:
            continue
        bond_objects.append((half, builder.add_bond(half, geometry)))

    snapshots = read_snapshots(settings.filepath, include_first=settings.include_first_snapshot)
    last_frame = 0
    if snapshots:
        last_frame = insert_atom_keyframes(atom_objects, snapshots, settings.slow_factor)
        last_frame = max(
            last_frame,
            insert_bond_keyframes(bond_objects, radii, snapshots, settings.slow_factor),
        )
        extend_scene_range(scene, first_frame=0, last_frame=last_frame)

    return LoadReport(
        filepath=structure.path,
        collection=builder.collection.name,
        atoms=len(atom_objects),
        bond_halves=len(bond_objects),
        snapshots=len(snapshots),
        last_frame=last_frame,
    )


def load(filepath, scene: bpy.types.Scene | None = None, **overrides: Any) -> LoadReport:
    """Load ``filepath``; ``overrides`` are :class:`LoaderSettings` fields."""
    return load_structure(LoaderSettings(filepath=Path(filepath), **overrides), scene=scene)


def load_from_environment(filepath=None, scene: bpy.types.Scene | None = None,
                          **overrides: Any) -> LoadReport:
    """Load the file named by the ``MOLECULAR_LOADER_*`` environment variables."""
    settings = LoaderSettings.from_environment(filepath=filepath, **overrides)
    return load_structure(settings, scene=scene)
