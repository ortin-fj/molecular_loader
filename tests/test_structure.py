"""Tests of the ASE reader, without Blender."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from molecular_loader.structure import read_snapshots, read_structure

DATA = Path(__file__).parent / "data"
ETHANOL = DATA / "ethanol_3frames.xyz"


def bonds(matrix: np.ndarray) -> int:
    """Number of bonded pairs, whichever direction the matrix lists them in."""
    matrix = np.asarray(matrix)
    symmetric = (matrix + matrix.T) != 0
    return int(np.triu(symmetric, k=1).sum())


def test_read_structure_reads_symbols_and_positions():
    structure = read_structure(ETHANOL)

    assert structure.number_of_atoms == 9
    assert structure.symbols == ("C", "C", "O", "H", "H", "H", "H", "H", "H")
    assert structure.positions.shape == (9, 3)


def test_read_structure_detects_the_bonds_of_ethanol():
    # ethanol: one C-C, one C-O, one O-H and five C-H bonds
    assert bonds(read_structure(ETHANOL).connectivity) == 8


def test_connectivity_lists_every_bond():
    connectivity = read_structure(ETHANOL).connectivity
    assert bonds(connectivity) == 8
    assert connectivity.shape == (9, 9)


def test_read_snapshots_keys_every_frame():
    snapshots = read_snapshots(ETHANOL)
    assert len(snapshots) == 3
    assert snapshots[0].shape == (9, 3)
    assert not np.allclose(snapshots[0], snapshots[-1])


def test_read_snapshots_can_skip_the_first_snapshot():
    snapshots = read_snapshots(ETHANOL, include_first=False)
    assert len(snapshots) == 2


def test_read_snapshots_of_a_single_geometry_is_empty():
    assert read_snapshots(DATA / "ethanol_single.xyz") == []


def test_edge_cases_file():
    """An element missing from the table, and two atoms at the same position."""
    structure = read_structure(DATA / "edge_cases.xyz")

    assert structure.symbols == ("Og", "H", "H")
    # Og-H and Og-H are real bonds; the two H sit on top of each other, which the
    # connectivity matrix reports as a bond too (the loader drops it later).
    assert bonds(structure.connectivity) == 3
    assert np.array_equal(structure.positions[1], structure.positions[2])
