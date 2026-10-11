"""Decide which pairs of atoms are bonded.  This module only needs numpy."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

from molecular_loader.elements import get_atom_properties
from molecular_loader.geometry import BondHalf, bond_radius


def find_bonds(connectivity, symbols: Sequence[str],
               explicit_pairs: Iterable[tuple[int, int]] = ()) -> list[BondHalf]:
    """Return the bond halves to draw, in creation order.

    Args:
        connectivity: Square matrix of bonds, such as the one built by ASE's
            ``NeighborList.get_connectivity_matrix``.  Only the upper triangle is
            used, so self-connections are ignored.
        symbols: Chemical symbol of every atom, in the same order as the matrix.
        explicit_pairs: Extra bonds to draw, as pairs of atom indices; they are
            added on top of the detected ones.

    Raises:
        ValueError: If the matrix does not match the symbols or an extra bond
            refers to an atom that does not exist.
    """
    matrix = np.asarray(connectivity).copy()
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"connectivity must be a square matrix, got shape {matrix.shape}")
    if matrix.shape[0] != len(symbols):
        raise ValueError(
            f"connectivity has {matrix.shape[0]} atoms but {len(symbols)} symbols were given"
        )

    for first, second in explicit_pairs:
        if not 0 <= first < len(symbols) or not 0 <= second < len(symbols):
            raise ValueError(
                f"extra bond ({first}, {second}) is out of range for {len(symbols)} atoms"
            )
        matrix[first, second] = matrix[second, first] = 1

    # The matrix may be symmetric or triangular (ASE returns one direction per
    # bond unless `bothways` is used), so both halves are made explicit.
    bond_matrix = matrix + matrix.T

    rows, cols = np.where(np.triu(bond_matrix, k=1) > 0)
    halves: list[BondHalf] = []
    for owner, other in zip(rows.tolist(), cols.tolist(), strict=True):
        same_element = symbols[owner] == symbols[other]
        radius = bond_radius(
            get_atom_properties(symbols[owner])["radius"],
            get_atom_properties(symbols[other])["radius"],
        )
        halves.append(
            BondHalf(f"bond{owner + 1}_{other + 1}", owner, other, radius, same_element)
        )
        if not same_element:
            halves.append(
                BondHalf(f"bond{other + 1}_{owner + 1}", other, owner, radius, same_element)
            )

    return halves
