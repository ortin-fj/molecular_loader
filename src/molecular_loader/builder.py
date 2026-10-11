"""Build the objects of a structure: spheres, bond cylinders and a collection.

All the objects of one element share a single mesh datablock, and all the
cylinders of one element and radius share another one, so a trajectory of a few
thousand atoms only needs a handful of meshes.  Cylinders are built one unit
long and stretched with their Z scale, which is also how the animation works.
"""

from __future__ import annotations

from collections.abc import Sequence

import bmesh
import bpy
import numpy as np

from molecular_loader.elements import get_atom_properties
from molecular_loader.geometry import BondHalf
from molecular_loader.materials import MaterialCache

SPHERE_SEGMENTS = 18
SPHERE_RINGS = 16
CYLINDER_SEGMENTS = 32
SPHERE_MESH_NAME = "sphere_{symbol}"
CYLINDER_MESH_NAME = "cylinder_{symbol}"
ATOM_NAME = "Atom{index}"
SUBDIVISION_MODIFIER = "Subdivision"


def shade_smooth(mesh: bpy.types.Mesh) -> None:
    """Mark every face of ``mesh`` as smooth, like Object > Shade Smooth."""
    for polygon in mesh.polygons:
        polygon.use_smooth = True


def create_collection(name: str) -> bpy.types.Collection:
    """Create a collection called ``name`` and link it to the scene.

    A new collection is always created; Blender adds a numeric suffix when the
    name is already taken, so an existing scene is never modified.
    """
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


class SceneBuilder:
    """Create the objects of one structure inside a new collection."""

    def __init__(self, collection_name: str, symbols: Sequence[str],
                 materials: MaterialCache) -> None:
        self.symbols = tuple(symbols)
        self.materials = materials
        self.collection = create_collection(collection_name)
        self._sphere_meshes: dict[str, bpy.types.Mesh] = {}
        self._cylinder_meshes: dict[tuple[str, float], bpy.types.Mesh] = {}

    def _link(self, obj: bpy.types.Object) -> bpy.types.Object:
        self.collection.objects.link(obj)
        return obj

    def sphere_mesh(self, symbol: str) -> bpy.types.Mesh:
        """Return the shared sphere mesh of ``symbol``."""
        if symbol not in self._sphere_meshes:
            mesh = bpy.data.meshes.new(SPHERE_MESH_NAME.format(symbol=symbol))
            bm = bmesh.new()
            bmesh.ops.create_uvsphere(
                bm,
                u_segments=SPHERE_SEGMENTS,
                v_segments=SPHERE_RINGS,
                radius=get_atom_properties(symbol)["radius"],
            )
            bm.to_mesh(mesh)
            bm.free()
            shade_smooth(mesh)
            mesh.materials.append(self.materials.get(symbol))
            self._sphere_meshes[symbol] = mesh
        return self._sphere_meshes[symbol]

    def cylinder_mesh(self, symbol: str, radius: float) -> bpy.types.Mesh:
        """Return the shared unit-length cylinder mesh of ``symbol`` and ``radius``."""
        key = (symbol, radius)
        if key not in self._cylinder_meshes:
            mesh = bpy.data.meshes.new(CYLINDER_MESH_NAME.format(symbol=symbol))
            bm = bmesh.new()
            bmesh.ops.create_cone(
                bm,
                cap_ends=True,
                cap_tris=False,
                segments=CYLINDER_SEGMENTS,
                radius1=radius,
                radius2=radius,
                depth=1.0,
            )
            bm.to_mesh(mesh)
            bm.free()
            shade_smooth(mesh)
            mesh.materials.append(self.materials.get(symbol))
            self._cylinder_meshes[key] = mesh
        return self._cylinder_meshes[key]

    def add_atom(self, index: int, position, symbol: str = "") -> bpy.types.Object:
        """Create the sphere of the atom ``index`` (0 based), centred on ``position``."""
        symbol = symbol or self.symbols[index]
        obj = bpy.data.objects.new(ATOM_NAME.format(index=index + 1), self.sphere_mesh(symbol))
        obj.location = tuple(np.asarray(position, dtype=float))
        obj.modifiers.new(name=SUBDIVISION_MODIFIER, type="SUBSURF")
        return self._link(obj)

    def add_bond(self, half: BondHalf, geometry) -> bpy.types.Object:
        """Create the cylinder of one bond half.

        The cylinder is stretched with its Z scale, so the radius of the mesh
        stays that of the element while the bond changes length over time.
        """
        centre, length, theta, phi = geometry
        symbol = self.symbols[half.owner]
        obj = bpy.data.objects.new(half.name, self.cylinder_mesh(symbol, half.radius))
        obj.location = tuple(np.asarray(centre, dtype=float))
        obj.rotation_euler[1] = theta
        obj.rotation_euler[2] = phi
        obj.scale[2] = length
        return self._link(obj)
