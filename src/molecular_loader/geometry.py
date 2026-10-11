"""Geometry of the bond cylinders.  This module only needs numpy."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

#: Shortest cylinder that is drawn. Atoms that overlap more than the sum of
#: their radii would give a negative length, so the length is clamped to this
#: value; atoms that sit closer than this are considered the same position.
MINIMUM_LENGTH = 1e-4


class BondHalf(NamedTuple):
    """One of the cylinders that a bond is drawn with.

    A bond between two different elements is drawn as two halves, each one
    coloured after the element it belongs to.  A bond between two equal elements
    is a single cylinder, whose ``owner`` is the lower index.
    """

    name: str
    owner: int
    other: int
    radius: float
    same_element: bool


def bond_radius(radius_a: float, radius_b: float) -> float:
    """Radius of the cylinders of a bond: half of the thinnest atom."""
    return min(radius_a, radius_b) / 2


def cylinder_axis(theta: float, phi: float) -> np.ndarray:
    """Return the direction a cylinder built along Z points to after ``theta``/``phi``.

    ``theta`` and ``phi`` are the Y and Z angles of an XYZ Euler rotation, which
    is what Blender applies to the bond objects.
    """
    return np.array([
        math.cos(phi) * math.sin(theta),
        math.sin(phi) * math.sin(theta),
        math.cos(theta),
    ])


def bond_geometry(positions: Sequence[Sequence[float]], owner: int, other: int,
                  radii: Sequence[float], same_element: bool,
                  minimum_length: float = MINIMUM_LENGTH):
    """Return ``(centre, length, theta, phi)`` for the cylinder of a bond half.

    The cylinder is built along the Z axis and then rotated with ``theta`` (Y)
    and ``phi`` (Z), so a unit-length cylinder scaled by ``length`` starts at
    ``positions[owner]`` and ends at ``positions[other]``.  ``None`` is returned
    when both atoms sit at the same position.
    """
    owner_position = np.asarray(positions[owner], dtype=float)
    delta = np.asarray(positions[other], dtype=float) - owner_position
    distance = float(np.linalg.norm(delta))
    if distance < minimum_length:
        return None

    if same_element:
        length = distance
    else:
        length = radii[owner] + (distance - radii[owner] - radii[other]) / 2
    length = max(length, minimum_length)

    centre = owner_position + delta * (length / (2 * distance))
    theta = math.acos(max(-1.0, min(1.0, float(delta[2] / distance))))
    phi = math.atan2(float(delta[1]), float(delta[0]))
    return centre, length, theta, phi
