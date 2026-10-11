"""Tests of the bond detection, without Blender."""

from __future__ import annotations

import numpy as np
import pytest

from molecular_loader.bonds import find_bonds
from molecular_loader.elements import ATOM_PROPERTIES, DEFAULT_RADIUS


def connectivity(size: int, *pairs: tuple[int, int], diagonal: bool = False) -> np.ndarray:
    """Build a connectivity matrix with ``pairs`` bonded (and optionally self-bonds)."""
    matrix = np.zeros((size, size), dtype=int)
    for first, second in pairs:
        matrix[first, second] = matrix[second, first] = 1
    if diagonal:
        np.fill_diagonal(matrix, 1)
    return matrix


def test_equal_elements_are_one_cylinder():
    halves = find_bonds(connectivity(2, (0, 1)), ("C", "C"))
    assert [half.name for half in halves] == ["bond1_2"]
    assert halves[0].same_element is True
    assert halves[0].owner == 0 and halves[0].other == 1


def test_different_elements_are_two_halves():
    halves = find_bonds(connectivity(2, (0, 1)), ("C", "O"))
    assert [half.name for half in halves] == ["bond1_2", "bond2_1"]
    assert [half.owner for half in halves] == [0, 1]
    assert all(half.same_element is False for half in halves)


def test_radius_is_half_of_the_thinnest_atom():
    halves = find_bonds(connectivity(2, (0, 1)), ("C", "H"))
    assert halves[0].radius == pytest.approx(ATOM_PROPERTIES["H"]["radius"] / 2)
    assert halves[1].radius == pytest.approx(ATOM_PROPERTIES["H"]["radius"] / 2)


def test_self_connections_are_ignored():
    halves = find_bonds(connectivity(2, (0, 1), diagonal=True), ("C", "C"))
    assert [half.name for half in halves] == ["bond1_2"]


def test_extra_bonds_are_added_on_top_of_the_detected_ones():
    halves = find_bonds(connectivity(3, (0, 1)), ("C", "C", "C"), explicit_pairs=[(1, 2)])
    assert sorted(half.name for half in halves) == ["bond1_2", "bond2_3"]


def test_extra_bond_of_different_elements_gives_two_halves():
    halves = find_bonds(connectivity(2), ("C", "O"), explicit_pairs=[(0, 1)])
    assert sorted(half.name for half in halves) == ["bond1_2", "bond2_1"]


def test_unknown_element_falls_back_to_the_default_radius():
    halves = find_bonds(connectivity(2, (0, 1)), ("Og", "Og"))
    assert halves[0].radius == pytest.approx(DEFAULT_RADIUS / 2)


def test_matrix_with_a_wrong_shape_is_rejected():
    with pytest.raises(ValueError, match="square"):
        find_bonds(np.zeros((2, 3)), ("C", "C"))


def test_matrix_with_a_wrong_size_is_rejected():
    with pytest.raises(ValueError, match="symbols"):
        find_bonds(np.zeros((3, 3)), ("C", "C"))


def test_extra_bond_out_of_range_is_rejected():
    with pytest.raises(ValueError, match="out of range"):
        find_bonds(connectivity(2), ("C", "C"), explicit_pairs=[(0, 5)])


def test_connectivity_can_be_triangular():
    """ASE fills one direction per bond; the lower triangle may be empty."""
    upper = connectivity(3, (0, 1), (1, 2))
    lower = np.tril(upper, k=-1).T
    for matrix in (upper, lower, upper + upper.T):
        halves = find_bonds(matrix, ("C", "C", "C"))
        assert sorted(half.name for half in halves) == ["bond1_2", "bond2_3"]


def test_connectivity_argument_is_not_modified():
    matrix = connectivity(2, (0, 1))
    find_bonds(matrix, ("C", "C"), explicit_pairs=[(0, 1)])
    assert matrix.tolist() == connectivity(2, (0, 1)).tolist()
