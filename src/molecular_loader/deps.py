"""Locate (and install) the third-party libraries that the loader needs.

Blender ships its own Python interpreter and does not add the user
``site-packages`` directory to ``sys.path``, so ``ase`` (and the ``scipy`` it
imports) are usually not importable from inside Blender.  Two directories are
supported for them:

* ``<add-on directory>/libs``, filled by the *Install dependencies* button of the
  add-on and by ``scripts/install_dependencies.py``;
* every directory listed in the ``MOLECULAR_LOADER_LIBS`` environment variable,
  separated by ``os.pathsep``.

The directories are *appended* to ``sys.path`` so that the packages Blender
ships (numpy, and therefore the ABI that scipy was built against) keep priority.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

#: Directory, next to this file, where the add-on installs its dependencies.
LIBRARY_DIRECTORY = "libs"

#: Environment variable with extra dependency directories (``os.pathsep`` separated).
LIBRARY_PATH_VARIABLE = "MOLECULAR_LOADER_LIBS"

#: Packages that have to be present inside Blender's Python.
RUNTIME_PACKAGES = ("ase", "scipy")


def library_directories() -> list[Path]:
    """Return the candidate directories for the runtime dependencies."""
    directories = [Path(__file__).resolve().parent / LIBRARY_DIRECTORY]
    extra = os.environ.get(LIBRARY_PATH_VARIABLE, "")
    directories.extend(Path(part) for part in extra.split(os.pathsep) if part)
    return directories


def add_library_directories() -> list[Path]:
    """Append the existing dependency directories to ``sys.path``.

    Idempotent, so it is safe to call it again after installing the packages.
    """
    added = []
    for directory in library_directories():
        if directory.is_dir() and str(directory) not in sys.path:
            sys.path.append(str(directory))
            added.append(directory)
    return added


def install_command(python, target, packages: Sequence[str] = RUNTIME_PACKAGES,
                    use_uv: bool = True) -> list[str]:
    """Return the command that installs ``packages`` into ``target``.

    ``--no-deps`` is deliberate: ``ase`` and ``scipy`` are installed without
    their own copy of numpy, so they use the (recent) numpy that Blender ships.
    """
    uv = shutil.which("uv") if use_uv else None
    if uv:
        return [uv, "pip", "install", "--python", str(python), "--target", str(target),
                "--no-deps", "--upgrade", *packages]
    return [str(python), "-m", "pip", "install", "--target", str(target),
            "--no-deps", "--upgrade", *packages]


def install_dependencies(target, python=None, packages: Iterable[str] = RUNTIME_PACKAGES,
                         use_uv: bool = True) -> subprocess.CompletedProcess:
    """Install ``packages`` into ``target`` and return the finished process.

    The command runs with the interpreter in ``python`` (Blender's own
    interpreter by default), so the wheels match the Blender version.
    """
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)
    python = Path(python) if python else Path(sys.executable)
    command = install_command(python, target, tuple(packages), use_uv=use_uv)
    return subprocess.run(command, capture_output=True, text=True)


def missing_packages(packages: Iterable[str] = RUNTIME_PACKAGES) -> list[str]:
    """Return the names of ``packages`` that cannot be imported."""
    import importlib.util

    add_library_directories()
    importlib.invalidate_caches()
    return [name for name in packages if importlib.util.find_spec(name) is None]
