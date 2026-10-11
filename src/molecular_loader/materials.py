"""One Cycles material per element.

The shader is a Principled BSDF, which covers Blender 3.x, 4.x and 5.x; the
glossy/diffuse mix of the original script is deprecated (``ShaderNodeBsdfGlossy``
and the removal of ``Material.use_nodes`` planned for Blender 6).
"""

from __future__ import annotations

import bpy

from molecular_loader.elements import get_atom_properties

#: Roughness of the element shaders: low enough to look glossy in Cycles.
ROUGHNESS = 0.25

MATERIAL_NAME = "{symbol}_mat"


def create_basic_material(mat_name: str, colour) -> bpy.types.Material:
    """Return the material called ``mat_name``, creating it if needed."""
    existing = bpy.data.materials.get(mat_name)
    if existing is not None:
        return existing

    material = bpy.data.materials.new(mat_name)
    material.diffuse_color = colour  # viewport colour, RGBA
    material.use_nodes = True

    nodes = material.node_tree.nodes
    nodes.clear()

    shader = nodes.new(type="ShaderNodeBsdfPrincipled")
    shader.inputs["Base Color"].default_value = colour
    shader.inputs["Roughness"].default_value = ROUGHNESS
    shader.location = -200, 0

    output = nodes.new(type="ShaderNodeOutputMaterial")
    output.location = 200, 0

    material.node_tree.links.new(shader.outputs["BSDF"], output.inputs["Surface"])
    return material


class MaterialCache:
    """The material of every element used by one structure."""

    def __init__(self) -> None:
        self._materials: dict[str, bpy.types.Material] = {}

    def get(self, symbol: str) -> bpy.types.Material:
        """Return the material of ``symbol``, creating it the first time."""
        if symbol not in self._materials:
            colour = get_atom_properties(symbol)["colour"]
            self._materials[symbol] = create_basic_material(
                MATERIAL_NAME.format(symbol=symbol), colour
            )
        return self._materials[symbol]

    @property
    def symbols(self) -> list[str]:
        return sorted(self._materials)
