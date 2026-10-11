"""Blender interface: File > Import, the sidebar panel, and their operators.

Everything the loader needs can be set here, so no file has to be edited before
loading a structure.
"""

from __future__ import annotations

from pathlib import Path

import bpy
from bpy.props import BoolProperty, IntProperty, PointerProperty, StringProperty
from bpy.types import Operator, Panel, PropertyGroup
from bpy_extras.io_utils import ImportHelper

from molecular_loader import deps
from molecular_loader.loader import load_structure
from molecular_loader.settings import LoaderSettings, parse_bond_pairs

PANEL_CATEGORY = "Molecular Loader"
STATUS_IDLE = "No structure loaded yet."
FILTER_GLOB = "*.xyz;*.pdb;*.cif;*.mol;*.sdf;*.traj;*.vasp;*.gen;*.cube;*.out;*.log"

_dependency_status = ""


def refresh_dependency_status() -> str:
    """Remember whether ASE can be imported, for the panel to display."""
    global _dependency_status
    missing = deps.missing_packages()
    if missing:
        _dependency_status = "Missing libraries: " + ", ".join(missing)
    else:
        import ase

        _dependency_status = f"ASE {ase.__version__} available"
    return _dependency_status


def settings_from_properties(properties: MolecularLoaderProperties) -> LoaderSettings:
    """Turn the panel properties into :class:`LoaderSettings`."""
    return LoaderSettings(
        filepath=Path(bpy.path.abspath(properties.filepath)),
        explicit_bonds=parse_bond_pairs(properties.explicit_bonds),
        slow_factor=properties.slow_factor,
        include_first_snapshot=properties.include_first_snapshot,
        collection_name=properties.collection_name,
    )


def load_and_report(operator: Operator, properties: MolecularLoaderProperties) -> set[str]:
    """Load the structure of ``properties`` and show the outcome of the operator."""
    try:
        report = load_structure(settings_from_properties(properties))
    except Exception as error:  # noqa: BLE001 - an add-on must not break the session
        properties.status = f"Error: {error}"
        operator.report({"ERROR"}, f"Molecular Loader: {error}")
        return {"CANCELLED"}

    properties.status = report.summary()
    operator.report({"INFO"}, report.summary())
    return {"FINISHED"}


class MolecularLoaderProperties(PropertyGroup):
    """Options of the panel, stored in the scene."""

    filepath: StringProperty(
        name="Structure file",
        description="File to load; any format that ASE can read (xyz, pdb, cif...)",
        subtype="FILE_PATH",
    )
    collection_name: StringProperty(
        name="Collection",
        description="Collection that will hold the objects (empty: name of the file)",
    )
    slow_factor: IntProperty(
        name="Slow factor",
        description="Insert a keyframe every N frames",
        default=1,
        min=1,
    )
    include_first_snapshot: BoolProperty(
        name="Keyframe first snapshot",
        description="Uncheck to use the first snapshot only as static geometry",
        default=True,
    )
    explicit_bonds: StringProperty(
        name="Extra bonds",
        description='Bonds that are not detected, as pairs of atom indices, e.g. "0-2, 0-3"',
    )
    status: StringProperty(name="Status", default=STATUS_IDLE)


class MOLECULAR_LOADER_OT_import(Operator, ImportHelper):
    """Load a molecule or a trajectory, choosing the file in the browser."""

    bl_idname = "molecular_loader.import_structure"
    bl_label = "Import Molecular Structure"
    bl_description = "Load spheres, bonds and keyframes from a structure file"
    bl_options = {"REGISTER", "UNDO"}

    filename_ext = ".xyz"
    filter_glob: StringProperty(default=FILTER_GLOB, options={"HIDDEN"})
    slow_factor: IntProperty(name="Slow factor", default=1, min=1)
    include_first_snapshot: BoolProperty(name="Keyframe first snapshot", default=True)
    explicit_bonds: StringProperty(name="Extra bonds")
    collection_name: StringProperty(name="Collection")

    def invoke(self, context, event):
        properties = context.scene.molecular_loader
        self.slow_factor = properties.slow_factor
        self.include_first_snapshot = properties.include_first_snapshot
        self.explicit_bonds = properties.explicit_bonds
        self.collection_name = properties.collection_name
        if properties.filepath:
            self.filepath = bpy.path.abspath(properties.filepath)
        return ImportHelper.invoke(self, context, event)

    def execute(self, context):
        properties = context.scene.molecular_loader
        properties.filepath = self.filepath
        properties.slow_factor = self.slow_factor
        properties.include_first_snapshot = self.include_first_snapshot
        properties.explicit_bonds = self.explicit_bonds
        properties.collection_name = self.collection_name
        return load_and_report(self, properties)


class MOLECULAR_LOADER_OT_load(Operator):
    """Load the file written in the panel."""

    bl_idname = "molecular_loader.load_structure"
    bl_label = "Load Structure"
    bl_description = "Load the structure file of the panel"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        properties = context.scene.molecular_loader
        if not properties.filepath:
            self.report({"ERROR"}, "Choose a structure file first")
            return {"CANCELLED"}
        return load_and_report(self, properties)


class MOLECULAR_LOADER_OT_install_dependencies(Operator):
    """Install ASE where Blender can import it (needs internet)."""

    bl_idname = "molecular_loader.install_dependencies"
    bl_label = "Install dependencies"
    bl_description = (
        f"Install {', '.join(deps.RUNTIME_PACKAGES)} in the add-on directory. "
        "Blender is busy while it downloads; the log goes to the system console"
    )

    def execute(self, context):
        target = deps.library_directories()[0]
        self.report({"INFO"}, f"Installing {', '.join(deps.RUNTIME_PACKAGES)} in {target}...")
        result = deps.install_dependencies(target)
        output = f"{result.stdout or ''}{result.stderr or ''}"
        refresh_dependency_status()

        if result.returncode != 0 or deps.missing_packages():
            print(output)
            last_line = output.strip().splitlines()[-1] if output.strip() else "unknown error"
            self.report({"ERROR"}, f"Installation failed: {last_line} (see the system console)")
            return {"CANCELLED"}

        self.report({"INFO"}, f"Installed {', '.join(deps.RUNTIME_PACKAGES)} in {target}")
        return {"FINISHED"}


class MOLECULAR_LOADER_OT_check_dependencies(Operator):
    """Check again whether ASE can be imported."""

    bl_idname = "molecular_loader.check_dependencies"
    bl_label = "Check libraries"
    bl_description = "Check whether ASE can be imported from Blender's Python"

    def execute(self, context):
        message = refresh_dependency_status()
        self.report({"INFO"}, message)
        return {"FINISHED"}


class MOLECULAR_LOADER_PT_panel(Panel):
    """Sidebar panel with the options of the loader."""

    bl_idname = "MOLECULAR_LOADER_PT_panel"
    bl_label = "Molecular Loader"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = PANEL_CATEGORY

    def draw(self, context):
        layout = self.layout
        properties = context.scene.molecular_loader

        layout.prop(properties, "filepath")
        row = layout.row(align=True)
        row.operator("molecular_loader.import_structure", text="Choose file...", icon="FILEBROWSER")
        row.operator("molecular_loader.load_structure", text="Load", icon="IMPORT")

        options = layout.box()
        options.label(text="Options")
        options.prop(properties, "slow_factor")
        options.prop(properties, "include_first_snapshot")
        options.prop(properties, "collection_name")
        options.prop(properties, "explicit_bonds")

        layout.label(text=properties.status, icon="INFO")

        libraries = layout.box()
        libraries.label(
            text=_dependency_status or refresh_dependency_status(),
            icon="CHECKMARK" if not deps.missing_packages() else "ERROR",
        )
        row = libraries.row(align=True)
        row.operator("molecular_loader.install_dependencies", icon="CONSOLE")
        row.operator("molecular_loader.check_dependencies", text="", icon="FILE_REFRESH")


CLASSES = (
    MolecularLoaderProperties,
    MOLECULAR_LOADER_OT_import,
    MOLECULAR_LOADER_OT_load,
    MOLECULAR_LOADER_OT_install_dependencies,
    MOLECULAR_LOADER_OT_check_dependencies,
    MOLECULAR_LOADER_PT_panel,
)


def menu_func_import(self, context):
    """Entry of File > Import."""
    self.layout.operator(
        MOLECULAR_LOADER_OT_import.bl_idname, text="Molecular Structure (ASE, xyz, pdb...)"
    )


def register() -> None:
    for cls in CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.Scene.molecular_loader = PointerProperty(type=MolecularLoaderProperties)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)
    refresh_dependency_status()


def unregister() -> None:
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)
    del bpy.types.Scene.molecular_loader
    for cls in reversed(CLASSES):
        bpy.utils.unregister_class(cls)
