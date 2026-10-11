"""Read structures and trajectories with ASE.

This is the only module that imports ASE, so the rest of the package can be
tested (and imported by Blender) without it.  ASE needs ``scipy`` as well, see
:mod:`molecular_loader.deps`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

_MISSING_ASE_MESSAGE = (
    "ASE is not available inside Blender's Python. Install the dependencies with "
    "the 'Install dependencies' button of the Molecular Loader panel, with "
    "scripts/install_dependencies.py, or manually with "
    "'uv pip install --python <blender python> --target <add-on>/libs --no-deps ase scipy'."
)

try:
    from ase.io import read
    from ase.neighborlist import NeighborList, natural_cutoffs
except ModuleNotFoundError as error:  # pragma: no cover - depends on the install
    raise ModuleNotFoundError(_MISSING_ASE_MESSAGE) from error


@dataclass(frozen=True)
class Structure:
    """The first snapshot of a file, ready to be built in the scene."""

    path: Path
    symbols: tuple[str, ...]
    positions: np.ndarray
    connectivity: np.ndarray

    @property
    def number_of_atoms(self) -> int:
        return len(self.symbols)


def read_structure(path) -> Structure:
    """Read ``path`` and return its atoms, positions and connectivity.

    Bonds are detected with ASE's covalent radii, which is what makes the
    ``explicit_bonds`` setting necessary for some metals.
    """
    atoms = read(path, index=0)
    neighbors = NeighborList(natural_cutoffs(atoms))
    neighbors.update(atoms)
    return Structure(
        path=Path(path),
        symbols=tuple(atoms.get_chemical_symbols()),
        positions=np.asarray(atoms.get_positions(), dtype=float),
        connectivity=neighbors.get_connectivity_matrix().toarray(),
    )


def read_snapshots(path, include_first: bool = True) -> list[np.ndarray]:
    """Return the positions of every snapshot that should be keyframed.

    A file with a single geometry returns an empty list, so a static structure
    gets no keyframes at all.  With ``include_first`` the first snapshot is
    keyframed too; otherwise it is left as the static geometry of the scene.
    """
    images = read(path, index=":")
    if len(images) <= 1:
        return []
    if not include_first:
        images = images[1:]
    return [np.asarray(image.get_positions(), dtype=float) for image in images]
