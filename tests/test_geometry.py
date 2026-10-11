"""Tests of the cylinder geometry, without Blender."""

from __future__ import annotations

import math

import numpy as np
import pytest

from molecular_loader.geometry import MINIMUM_LENGTH, bond_geometry, bond_radius, cylinder_axis

RADII = np.array([0.5, 0.5, 0.25])


def direction(*args) -> np.ndarray:
    """Unit vector of the first two components of ``args``."""
    vector = np.asarray(args, dtype=float)
    return vector / np.linalg.norm(vector)


@pytest.mark.parametrize(
    "positions",
    [
        [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        [(1.0, 1.0, 1.0), (1.0, 3.0, 1.0)],
        [(0.0, 0.0, 0.0), (0.0, 0.0, 2.0)],
        [(0.0, 0.0, 0.0), (-1.5, 2.0, -0.75)],
    ],
)
def test_cylinder_axis_matches_the_bond(positions):
    centre, length, theta, phi = bond_geometry(positions, 0, 1, RADII, same_element=True)
    expected = direction(*(np.asarray(positions[1], dtype=float) - positions[0]))
    assert np.allclose(cylinder_axis(theta, phi), expected, atol=1e-12)
    assert centre == pytest.approx(np.asarray(positions[0]) + expected * length / 2)


@pytest.mark.parametrize("positions_length", [1.0, 2.5])
def test_equal_elements_span_the_two_atoms(positions_length):
    positions = [(0.0, 0.0, 0.0), (positions_length, 0.0, 0.0)]
    centre, length, theta, phi = bond_geometry(positions, 0, 1, RADII, same_element=True)
    assert length == pytest.approx(positions_length)
    assert centre == pytest.approx(np.array([positions_length / 2, 0.0, 0.0]))
    assert theta == pytest.approx(math.pi / 2)
    assert phi == pytest.approx(0.0)


def test_the_two_halves_of_a_bond_meet_in_the_middle():
    positions = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    radii = np.array([0.4, 0.1])
    first = bond_geometry(positions, 0, 1, radii, same_element=False)
    second = bond_geometry(positions, 1, 0, radii, same_element=False)
    assert first[1] + second[1] == pytest.approx(2.0)


def test_a_half_only_reaches_half_way_through_the_gap():
    positions = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    radii = np.array([0.4, 0.1])
    centre, length, _, _ = bond_geometry(positions, 0, 1, radii, same_element=False)
    # radius of the first atom plus half of the distance that is not covered by
    # either atom
    assert length == pytest.approx(0.4 + (2.0 - 0.4 - 0.1) / 2)
    assert centre[0] == pytest.approx(length / 2)


def test_coincident_atoms_give_no_cylinder():
    positions = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
    assert bond_geometry(positions, 0, 1, RADII, same_element=True) is None


def test_negative_lengths_are_clamped():
    # The second atom is thinner, so the formula gives a negative length
    positions = [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)]
    radii = np.array([0.1, 1.0])
    centre, length, _, _ = bond_geometry(positions, 0, 1, radii, same_element=False)
    assert length == pytest.approx(MINIMUM_LENGTH)
    assert centre[0] == pytest.approx(MINIMUM_LENGTH / 2)


def test_bond_radius_is_half_of_the_thinnest_element():
    assert bond_radius(0.5, 0.2) == pytest.approx(0.1)
    assert bond_radius(0.2, 0.5) == pytest.approx(0.1)


def test_cylinder_axis_is_a_unit_vector():
    for theta, phi in [(0.0, 0.0), (math.pi / 2, math.pi / 3), (2.5, -1.2)]:
        assert np.linalg.norm(cylinder_axis(theta, phi)) == pytest.approx(1.0)
