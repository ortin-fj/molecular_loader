"""Load a structure into Blender from the command line.

    blender --background --python scripts/load_in_blender.py -- molecule.xyz
    blender --python scripts/load_in_blender.py -- molecule.xyz --save-blend scene.blend

Options left out are taken from the ``MOLECULAR_LOADER_*`` environment variables,
so an existing pipeline can drive the loader without touching any file:

    MOLECULAR_LOADER_FILE=molecule.xyz MOLECULAR_LOADER_SLOW_FACTOR=5 \\
        blender --background --python scripts/load_in_blender.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from molecular_loader.deps import LIBRARY_PATH_VARIABLE, add_library_directories  # noqa: E402
from molecular_loader.settings import (  # noqa: E402
    parse_bond_pairs,
)

#: Flags of Blender's own command line, used to detect that the script is being
#: run by Blender without a ``--`` separator.
BLENDER_FLAGS = ("-b", "--background", "-P", "--python", "--python-text",
                 "--python-expr", "--python-console")


def script_arguments(argv: list[str]) -> list[str]:
    """Return the arguments after the ``--`` that separates them from Blender's.

    Blender passes its own flags to the interpreter, so everything up to ``--``
    is ignored there; when the script is run by a plain interpreter those are
    the only arguments there are.
    """
    if "--" in argv:
        return argv[argv.index("--") + 1:]
    if any(flag in argv for flag in BLENDER_FLAGS):
        return []
    return argv[1:]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("file", nargs="?", help="structure file; defaults to "
                                               "MOLECULAR_LOADER_FILE")
    parser.add_argument("--slow-factor", type=int, default=None,
                        help="insert a keyframe every N frames")
    parser.add_argument("--explicit-bonds", default=None,
                        help='bonds that are not detected, e.g. "0-2, 0-3"')
    parser.add_argument("--collection", default=None,
                        help="collection that holds the new objects")
    parser.add_argument("--skip-first-snapshot", action="store_true",
                        help="use the first snapshot only as static geometry")
    parser.add_argument("--deps", type=Path, action="append", default=[],
                        help="directory with ase/scipy, like " + LIBRARY_PATH_VARIABLE)
    parser.add_argument("--save-blend", type=Path, default=None,
                        help="save the scene to this .blend file")
    return parser


def main(argv=None) -> int:
    arguments = build_parser().parse_args(script_arguments(argv if argv is not None else sys.argv))

    for directory in arguments.deps:  # has to happen before ASE is imported
        sys.path.append(str(directory.expanduser().resolve()))
    add_library_directories()

    overrides = {}
    if arguments.slow_factor is not None:
        overrides["slow_factor"] = arguments.slow_factor
    if arguments.explicit_bonds is not None:
        overrides["explicit_bonds"] = parse_bond_pairs(arguments.explicit_bonds)
    if arguments.collection is not None:
        overrides["collection_name"] = arguments.collection
    if arguments.skip_first_snapshot:
        overrides["include_first_snapshot"] = False

    try:
        # Imported here because it needs ASE (and Blender itself, so this script
        # can only run inside Blender).
        from molecular_loader.loader import load_from_environment

        report = load_from_environment(filepath=arguments.file, **overrides)
    except ModuleNotFoundError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        print("Tip: inside Blender, write the options after the '--' separator, "
              "for example: blender --background --python scripts/load_in_blender.py "
              "-- molecule.xyz", file=sys.stderr)
        return 1

    print(report.summary())

    if arguments.save_blend:
        import bpy

        path = arguments.save_blend.expanduser().resolve()
        bpy.ops.wm.save_as_mainfile(filepath=str(path))
        print(f"Saved {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
