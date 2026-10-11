"""Install ASE where Blender can import it.

Blender ships its own Python interpreter, and it is not the interpreter of the
project's virtual environment, so this installs ASE and scipy *for Blender*
(with ``--no-deps``, so they reuse the numpy that Blender already ships).

Examples:
    uv run python scripts/install_dependencies.py            # into src/molecular_loader/libs
    python scripts/install_dependencies.py --blender /usr/bin/blender
    python scripts/install_dependencies.py --python /path/to/python --target /tmp/libs
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from molecular_loader.deps import RUNTIME_PACKAGES, install_dependencies  # noqa: E402

DEFAULT_TARGET = REPOSITORY / "src" / "molecular_loader" / "libs"


def blender_python(blender: str) -> Path:
    """Return the interpreter that ``blender`` executes (its own bundled Python)."""
    try:
        result = subprocess.run(
            [blender, "--background", "--factory-startup", "--python-expr",
             "import sys; print(sys.executable)"],
            capture_output=True, text=True, check=False,
        )
    except FileNotFoundError as error:
        raise SystemExit(f"Blender not found: {error}") from error

    candidates = [
        Path(line.strip())
        for line in result.stdout.splitlines()
        if line.strip().startswith("/") and "python" in Path(line.strip()).name.lower()
    ]
    if result.returncode != 0 or not candidates:
        raise SystemExit(
            f"could not find the Python interpreter of {blender}:\n{result.stdout}\n{result.stderr}"
        )
    return candidates[-1]


def verify(python: Path, target: Path) -> subprocess.CompletedProcess:
    """Import the packages from ``target`` with ``python``, to be sure they work."""
    return subprocess.run(
        [str(python), "-c",
         f"import sys; sys.path.append({str(target)!r}); import ase, scipy; "
         "print('ase', ase.__version__, '| scipy', scipy.__version__)"],
        capture_output=True, text=True, check=False,
    )


def main(arguments=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET,
                        help=f"directory to install into (default: {DEFAULT_TARGET})")
    interpreter = parser.add_mutually_exclusive_group()
    interpreter.add_argument("--python", type=Path,
                             help="interpreter to install for (default: Blender's own)")
    interpreter.add_argument("--blender", default=shutil.which("blender") or "blender",
                             help="Blender executable used to find its Python")
    parser.add_argument("--no-uv", action="store_true", help="use 'python -m pip' instead of uv")
    options = parser.parse_args(arguments)

    python = options.python or blender_python(options.blender)
    target = options.target.expanduser().resolve()
    print(f"Installing {', '.join(RUNTIME_PACKAGES)} into {target}")
    print(f"Using {python}")

    result = install_dependencies(target, python=python, use_uv=not options.no_uv)
    print(result.stdout.strip() or result.stderr.strip())
    if result.returncode != 0:
        return result.returncode

    check = verify(python, target)
    print(check.stdout.strip() or check.stderr.strip())
    if check.returncode != 0:
        print("The packages were installed but Blender could not import them.", file=sys.stderr)
    else:
        print(f"Done. Add {target} to the MOLECULAR_LOADER_LIBS variable, or install the "
              "add-on from that directory, to make it importable from Blender.")
    return check.returncode


if __name__ == "__main__":
    raise SystemExit(main())
