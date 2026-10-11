"""Build the zip that installs Molecular Loader in Blender.

    uv run python scripts/build_addon.py

The archive holds the package together with ``blender_manifest.toml``, so the
same zip is valid as an extension (Blender 4.2 and newer) and as a legacy add-on
(Blender 3.6 - 4.1).  Install it from *Edit > Preferences > Add-ons > Install
from Disk*.
"""

from __future__ import annotations

import shutil
import subprocess
import zipfile
from collections.abc import Iterator
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parent.parent
PACKAGE = REPOSITORY / "src" / "molecular_loader"
MANIFEST = REPOSITORY / "blender_manifest.toml"
DIST = REPOSITORY / "dist"

#: Not part of the archive: installed dependencies and Python caches.
EXCLUDED_DIRECTORIES = {"libs", "__pycache__"}


def read_version() -> str:
    """Read the version from the manifest, which is the single source of truth."""
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.startswith("version"):
            return line.split("=", 1)[1].strip().strip('"')
    raise SystemExit(f"no version found in {MANIFEST}")


def package_files() -> Iterator[tuple[Path, Path]]:
    """Yield ``(source, name inside the archive)`` for every file to bundle."""
    for path in sorted(PACKAGE.rglob("*")):
        relative = path.relative_to(PACKAGE)
        if EXCLUDED_DIRECTORIES.intersection(relative.parts) or path.suffix == ".pyc":
            continue
        if path.is_file():
            yield path, Path(PACKAGE.name) / relative


def build(validate: bool = True) -> Path:
    version = read_version()
    DIST.mkdir(exist_ok=True)
    archive = DIST / f"{PACKAGE.name}-{version}.zip"

    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
        bundle.write(MANIFEST, Path(PACKAGE.name) / MANIFEST.name)
        for source, name in package_files():
            bundle.write(source, name)

    print(f"{archive} ({archive.stat().st_size / 1024:.0f} kB)")

    blender = shutil.which("blender")
    if validate and blender:
        result = subprocess.run([blender, "--command", "extension", "validate", str(archive)],
                                capture_output=True, text=True, check=False)
        print(result.stdout.strip() or result.stderr.strip())
        if result.returncode != 0:
            raise SystemExit("Blender rejected the archive")
    elif validate:
        print("Blender was not found, skipping the validation")
    return archive


if __name__ == "__main__":
    build()
