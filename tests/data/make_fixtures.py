"""Regenerate the structure files used by the tests.

    uv run python tests/data/make_fixtures.py

``ethanol_3frames.xyz`` is ethanol (C-C, C-O, O-H and five C-H bonds) rotating
and expanding slightly, so both the single-cylinder (equal elements) and the
two-halves (different elements) cases are covered.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase.build import molecule
from ase.io import write

DATA = Path(__file__).parent


def rotation_z(angle_degrees: float) -> np.ndarray:
    angle = np.radians(angle_degrees)
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]])


def main() -> None:
    base = molecule("CH3CH2OH")
    base.center(vacuum=2.0)

    frames = []
    for index in range(3):
        image = base.copy()
        image.positions = base.positions @ rotation_z(15.0 * index).T * (1.0 + 0.03 * index)
        frames.append(image)

    write(DATA / "ethanol_3frames.xyz", frames)
    write(DATA / "ethanol_single.xyz", [frames[0]])
    print(f"wrote ethanol_3frames.xyz and ethanol_single.xyz ({len(frames)} frames)")


if __name__ == "__main__":
    main()
