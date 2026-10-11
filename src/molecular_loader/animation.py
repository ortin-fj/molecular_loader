"""Write the animation of a trajectory.

Cylinders are one unit long, so their Z scale is the length of the bond, which
makes the scale keyframes independent of the geometry they were created with.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import bpy
import numpy as np

from molecular_loader.geometry import BondHalf, bond_geometry


def insert_atom_keyframes(atom_objects: Sequence[bpy.types.Object],
                          snapshots: Iterable[np.ndarray], slow_factor: int = 1) -> int:
    """Keyframe the position of every atom and return the last frame used."""
    last_frame = 0
    for index, positions in enumerate(snapshots):
        frame = index * slow_factor
        # A trajectory with fewer atoms in a later snapshot keeps the previous
        # pose of the missing ones instead of failing the whole load.
        for obj, position in zip(atom_objects, positions, strict=False):
            obj.location = tuple(position)
            obj.keyframe_insert(data_path="location", frame=frame)
        last_frame = frame
    return last_frame


def insert_bond_keyframes(bond_objects: Sequence[tuple[BondHalf, bpy.types.Object]],
                          radii: Sequence[float], snapshots: Iterable[np.ndarray],
                          slow_factor: int = 1) -> int:
    """Keyframe location, rotation and scale of every bond half; returns the last frame."""
    last_frame = 0
    for index, positions in enumerate(snapshots):
        frame = index * slow_factor
        for half, obj in bond_objects:
            if half.owner >= len(positions) or half.other >= len(positions):
                continue
            geometry = bond_geometry(positions, half.owner, half.other, radii, half.same_element)
            if geometry is None:
                continue
            centre, length, theta, phi = geometry
            obj.location = tuple(centre)
            obj.rotation_euler[1] = theta
            obj.rotation_euler[2] = phi
            obj.scale[2] = length
            for path in ("location", "rotation_euler", "scale"):
                obj.keyframe_insert(data_path=path, frame=frame)
        last_frame = frame
    return last_frame


def extend_scene_range(scene: bpy.types.Scene, first_frame: int, last_frame: int) -> None:
    """Grow the frame range of the scene so the new keyframes are visible."""
    scene.frame_start = min(scene.frame_start, first_frame)
    scene.frame_end = max(scene.frame_end, last_frame)
