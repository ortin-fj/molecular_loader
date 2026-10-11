"""End-to-end test of the loader, to be run *inside* Blender.

    blender --background --factory-startup --python tests/blender/smoke_test.py

It builds the scene of the fixtures in ``tests/data`` and checks the result
against expectations that are computed here: object names, materials, shared
meshes, keyframes and the world position of every cylinder.  The add-on is
enabled as well, to test the File > Import operator.

ASE (and scipy) must be importable from Blender's Python; the script looks for
them in ``src/molecular_loader/libs`` and in ``MOLECULAR_LOADER_LIBS``.
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "src"))
LIBRARIES = REPOSITORY / "src" / "molecular_loader" / "libs"
if LIBRARIES.is_dir():
    sys.path.append(str(LIBRARIES))

import bpy  # noqa: E402
import numpy as np  # noqa: E402
from mathutils import Vector  # noqa: E402

DATA = REPOSITORY / "tests" / "data"
TRAJECTORY = DATA / "ethanol_3frames.xyz"
SINGLE = DATA / "ethanol_single.xyz"
EDGE_CASES = DATA / "edge_cases.xyz"

FAILURES: list[str] = []
CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        FAILURES.append(message)


def fresh_scene():
    """An empty scene, without the default cube, camera and light."""
    bpy.ops.wm.read_factory_settings(use_empty=True)
    return bpy.context.scene


def get_fcurves(obj):
    """F-curves of an object, across the legacy and the slotted action APIs."""
    action = obj.animation_data.action if obj.animation_data else None
    if action is None:
        return []
    if hasattr(action, "fcurves"):
        return list(action.fcurves)
    curves = []
    for layer in getattr(action, "layers", []):
        for strip in getattr(layer, "strips", []):
            for channelbag in getattr(strip, "channelbags", []):
                curves.extend(channelbag.fcurves)
    return curves


def keyframe_count(obj) -> int:
    return sum(len(curve.keyframe_points) for curve in get_fcurves(obj))


def operator_available(module: str, name: str) -> bool:
    """Whether ``bpy.ops.<module>.<name>`` is registered.

    ``hasattr`` cannot be used: the ``bpy.ops`` namespace answers any name and
    only fails when the operator is actually called.
    """
    try:
        getattr(getattr(bpy.ops, module), name).get_rna_type()
    except (AttributeError, KeyError):
        return False
    return True


def mesh_depth(mesh) -> float:
    """Extent of a mesh along its local Z axis."""
    coordinates = [vertex.co.z for vertex in mesh.vertices]
    return max(coordinates) - min(coordinates)


def endpoints(obj) -> list[np.ndarray]:
    """World position of both ends of a (one unit long) bond cylinder."""
    return [
        np.array(obj.matrix_world @ Vector((0.0, 0.0, sign * 0.5))) for sign in (-1.0, 1.0)
    ]


def test_trajectory() -> None:
    """A three-snapshot trajectory, keyframed from the first snapshot."""
    from ase.io import read

    from molecular_loader import LoaderSettings, load_structure

    scene = fresh_scene()
    report = load_structure(LoaderSettings(filepath=TRAJECTORY), scene=scene)

    frames = read(TRAJECTORY, index=":")
    symbols = frames[0].get_chemical_symbols()

    check(report.atoms == 9, f"expected 9 atoms, got {report.atoms}")
    check(report.collection == "ethanol_3frames", f"collection: {report.collection}")
    check(report.snapshots == 3, f"expected 3 keyframed snapshots, got {report.snapshots}")
    check(report.last_frame == 2, f"last frame: {report.last_frame}")

    collection = bpy.data.collections[report.collection]
    names = {obj.name for obj in collection.objects}
    check(len(names) == len(collection.objects), "object names are not unique")

    expected_bonds = {
        f"bond{first + 1}_{second + 1}"
        for first in range(len(symbols))
        for second in range(len(symbols))
        if first != second
        and (symbols[first] != symbols[second] or first < second)
        and _bonded(frames[0], first, second)
    }
    expected_names = {f"Atom{index + 1}" for index in range(9)} | expected_bonds
    check(
        names == expected_names,
        f"unexpected object names: {sorted(names.symmetric_difference(expected_names))[:6]}",
    )
    check(len(expected_bonds) == 15, f"fixture changed: {len(expected_bonds)} bond halves")

    # materials: one per element, and every object uses that of its element
    for index, symbol in enumerate(symbols):
        atom = bpy.data.objects[f"Atom{index + 1}"]
        materials = [material.name for material in atom.data.materials]
        check(materials == [f"{symbol}_mat"], f"Atom{index + 1} material: {materials}")
    for name in expected_bonds:
        owner = int(name.split("_")[0][4:]) - 1
        materials = [material.name for material in bpy.data.objects[name].data.materials]
        check(materials == [f"{symbols[owner]}_mat"], f"{name} material: {materials}")

    # meshes are shared: one sphere per element, one cylinder per element and radius
    spheres = {mesh.name for mesh in bpy.data.meshes if mesh.name.startswith("sphere_")}
    check(len(spheres) == 3, f"expected 3 shared sphere meshes, got {sorted(spheres)}")
    cylinders = {mesh.name: mesh for mesh in bpy.data.meshes if mesh.name.startswith("cylinder_")}
    check(
        len(cylinders) == 5,
        f"expected 5 shared cylinder meshes (element and radius), got {sorted(cylinders)}",
    )
    check(
        all(mesh_depth(mesh) == 1.0 for mesh in cylinders.values()),
        "cylinder meshes must be one unit long, the objects scale them",
    )

    # the scene: subdivision only on the atoms, everything smooth, everything animated
    check(
        {obj.name for obj in collection.objects if obj.modifiers.get("Subdivision")}
        == {f"Atom{index + 1}" for index in range(9)},
        "subdivision surface should only be added to the atoms",
    )
    check(
        all(
            all(polygon.use_smooth for polygon in obj.data.polygons)
            for obj in collection.objects
        ),
        "every face of every object should be smooth",
    )
    for name in names:
        obj = bpy.data.objects[name]
        channels = 3 if name.startswith("Atom") else 9
        check(
            keyframe_count(obj) == channels * 3,
            f"{name} has {keyframe_count(obj)} keyframes, expected {channels * 3}",
        )
    check(scene.frame_start == 0, f"scene frame start: {scene.frame_start}")
    check(scene.frame_end >= 2, f"scene frame end: {scene.frame_end}")

    # geometry: at every keyframed frame the bonds must touch their atoms
    radii = np.array([_radius(symbol) for symbol in symbols])
    for frame_index, image in enumerate(frames):
        scene.frame_set(frame_index)
        positions = image.get_positions()
        for index in range(9):
            location = np.array(bpy.data.objects[f"Atom{index + 1}"].matrix_world.translation)
            check(
                np.allclose(location, positions[index], atol=1e-5),
                f"Atom{index + 1} is misplaced at frame {frame_index}",
            )
        for name in expected_bonds:
            owner = int(name.split("_")[0][4:]) - 1
            other = int(name.split("_")[1]) - 1
            start, end = endpoints(bpy.data.objects[name])
            owner_position = positions[owner]
            other_position = positions[other]

            distance_to_owner = min(np.linalg.norm(start - owner_position),
                                    np.linalg.norm(end - owner_position))
            check(
                distance_to_owner < 1e-5,
                f"{name} does not start at Atom{owner + 1} at frame {frame_index}",
            )

            far = end if np.linalg.norm(start - owner_position) < 1e-5 else start
            direction = other_position - owner_position
            length = float(np.linalg.norm(direction))
            unit = direction / length
            along = float(np.dot(far - owner_position, unit))
            across = float(np.linalg.norm((far - owner_position) - along * unit))
            expected_length = (
                length if symbols[owner] == symbols[other]
                else radii[owner] + (length - radii[owner] - radii[other]) / 2
            )
            check(across < 1e-5, f"{name} is off axis at frame {frame_index} ({across})")
            check(
                abs(along - expected_length) < 1e-5,
                f"{name} is {along} long instead of {expected_length} at frame {frame_index}",
            )

    # the first snapshot is keyframed too, so frame 0 is the first snapshot
    scene.frame_set(0)
    check(
        np.allclose(
            np.array(bpy.data.objects["Atom1"].matrix_world.translation),
            frames[0].get_positions()[0],
            atol=1e-5,
        ),
        "frame 0 should hold the first snapshot",
    )


def test_first_snapshot_can_be_skipped() -> None:
    """The old behaviour: the first snapshot is only static geometry."""
    from molecular_loader import LoaderSettings, load_structure

    scene = fresh_scene()
    report = load_structure(
        LoaderSettings(filepath=TRAJECTORY, include_first_snapshot=False), scene=scene
    )
    scene.frame_set(0)

    from ase.io import read

    image = read(TRAJECTORY, index=":")[1]
    check(report.snapshots == 2, f"expected 2 snapshots, got {report.snapshots}")
    check(
        np.allclose(
            np.array(bpy.data.objects["Atom1"].matrix_world.translation),
            image.get_positions()[0],
            atol=1e-5,
        ),
        "frame 0 should hold the second snapshot when the first one is skipped",
    )


def test_single_geometry_has_no_keyframes() -> None:
    from molecular_loader import LoaderSettings, load_structure

    scene = fresh_scene()
    report = load_structure(LoaderSettings(filepath=SINGLE), scene=scene)

    check(report.snapshots == 0, f"expected no snapshots, got {report.snapshots}")
    check(
        all(obj.animation_data is None for obj in bpy.data.objects),
        "a single geometry must not be animated",
    )
    check(len(bpy.data.objects) == 24, f"objects: {len(bpy.data.objects)}")


def test_slow_factor_spreads_the_snapshots() -> None:
    from molecular_loader import LoaderSettings, load_structure

    scene = fresh_scene()
    report = load_structure(LoaderSettings(filepath=TRAJECTORY, slow_factor=10), scene=scene)

    check(report.last_frame == 20, f"last frame: {report.last_frame}")
    check(
        keyframe_count(bpy.data.objects["Atom1"]) == 3 * 3,
        "one keyframe per snapshot, in the three location channels",
    )


def test_extra_bonds_and_edge_cases() -> None:
    from molecular_loader import LoaderSettings, load_structure

    scene = fresh_scene()
    report = load_structure(
        LoaderSettings(filepath=TRAJECTORY, explicit_bonds=((5, 8), (4, 7))), scene=scene
    )
    check(report.bond_halves == 17, f"extra bonds: {report.bond_halves} bond halves")

    scene = fresh_scene()
    report = load_structure(LoaderSettings(filepath=EDGE_CASES), scene=scene)
    check(report.atoms == 3, f"edge case atoms: {report.atoms}")
    check(report.bond_halves == 4, f"edge case bond halves: {report.bond_halves}")
    names = {obj.name for obj in bpy.data.collections[report.collection].objects}
    check(
        not {"bond2_3", "bond3_2"} & names,
        f"atoms at the same position must not be bonded: {sorted(names)}",
    )
    check("Og_mat" in bpy.data.materials, "the unknown element should get a fallback material")

    scene = fresh_scene()
    report = load_structure(
        LoaderSettings(filepath=EDGE_CASES, collection_name="custom"), scene=scene
    )
    check(report.collection == "custom", f"collection: {report.collection}")


def test_addon_operators() -> None:
    """The add-on itself: registration and the File > Import operator."""
    import addon_utils

    # Loading factory settings disables every add-on, so it happens first.
    scene = fresh_scene()
    check(
        not operator_available("molecular_loader", "import_structure"),
        "the operator should not exist before the add-on is enabled",
    )

    module = addon_utils.enable("molecular_loader", default_set=True)
    check(module is not None, "the add-on could not be enabled")
    check(
        operator_available("molecular_loader", "import_structure"),
        "the import operator is missing after enabling the add-on",
    )
    check(hasattr(bpy.types.Scene, "molecular_loader"), "panel properties missing")

    result = bpy.ops.molecular_loader.import_structure(filepath=str(TRAJECTORY))
    check(result == {"FINISHED"}, f"operator returned {result}")
    check(len(bpy.data.collections) == 1, "the operator should create one collection")
    check(
        scene.molecular_loader.status.startswith("Loaded 9 atoms"),
        f"panel status: {scene.molecular_loader.status}",
    )
    check(
        scene.molecular_loader.filepath == str(TRAJECTORY),
        "the operator should remember the file",
    )

    # From Python an operator that reports an error raises instead of returning
    # CANCELLED, which is what the panel shows to the user.
    try:
        bpy.ops.molecular_loader.import_structure(filepath=str(DATA / "does-not-exist.xyz"))
        check(False, "loading a missing file should fail")
    except RuntimeError as error:
        check("No such file" in str(error), f"unexpected error: {error}")
    check(
        scene.molecular_loader.status.startswith("Error:"),
        f"panel status: {scene.molecular_loader.status}",
    )

    # Without a file the operator reports the problem instead of loading.
    # (Loading factory settings here would disable the add-on, so the scene
    # properties are cleared instead.)
    scene.molecular_loader.filepath = ""
    try:
        bpy.ops.molecular_loader.load_structure()
        check(False, "loading without a file should fail")
    except RuntimeError as error:
        check("Choose a structure file" in str(error), f"unexpected error: {error}")

    addon_utils.disable("molecular_loader")
    check(
        not operator_available("molecular_loader", "import_structure"),
        "the operator should be gone after disabling the add-on",
    )


def _bonded(image, first: int, second: int) -> bool:
    """Whether two atoms are bonded; ASE fills only one direction of the matrix."""
    from ase.neighborlist import NeighborList, natural_cutoffs

    neighbors = NeighborList(natural_cutoffs(image))
    neighbors.update(image)
    connectivity = neighbors.get_connectivity_matrix().toarray()
    return bool(connectivity[first, second] or connectivity[second, first])


def _radius(symbol: str) -> float:
    from molecular_loader.elements import get_atom_properties

    return get_atom_properties(symbol)["radius"]


def main() -> int:
    from molecular_loader import deps

    if deps.missing_packages():
        print(f"SKIP: install {', '.join(deps.RUNTIME_PACKAGES)} first "
              f"(missing {', '.join(deps.missing_packages())}). "
              "See the README, or use the 'Install dependencies' button.")
        return 2

    for test in (
        test_trajectory,
        test_first_snapshot_can_be_skipped,
        test_single_geometry_has_no_keyframes,
        test_slow_factor_spreads_the_snapshots,
        test_extra_bonds_and_edge_cases,
        test_addon_operators,
    ):
        try:
            test()
            print(f"ran {test.__name__}")
        except Exception:  # noqa: BLE001 - report every failure, keep testing
            FAILURES.append(f"{test.__name__} raised an exception")
            traceback.print_exc()

    print(f"\n{CHECKS} checks, {len(FAILURES)} failures")
    for failure in FAILURES:
        print(f"  FAILED: {failure}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
