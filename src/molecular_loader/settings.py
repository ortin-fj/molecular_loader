"""Settings that describe one load into Blender.

The same options can be given to the API (:class:`LoaderSettings`), typed in the
add-on panel, passed on the command line (``scripts/load_in_blender.py``) or set
as environment variables, so a scene can be produced without editing any file.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

#: Matches two integers separated by anything that is not a digit, which covers
#: "0-2", "0, 2", "0 2", "[[0, 2], [0, 3]]" and similar notations.
BOND_PAIR_PATTERN = re.compile(r"(\d+)\D+(\d+)")

FILE_VARIABLE = "MOLECULAR_LOADER_FILE"
SLOW_FACTOR_VARIABLE = "MOLECULAR_LOADER_SLOW_FACTOR"
EXPLICIT_BONDS_VARIABLE = "MOLECULAR_LOADER_EXPLICIT_BONDS"
COLLECTION_VARIABLE = "MOLECULAR_LOADER_COLLECTION"
INCLUDE_FIRST_VARIABLE = "MOLECULAR_LOADER_INCLUDE_FIRST_SNAPSHOT"

_TRUE = {"1", "true", "yes", "y", "on", "si", "sí"}
_FALSE = {"0", "false", "no", "n", "off"}


def parse_bond_pairs(text) -> tuple[tuple[int, int], ...]:
    """Read extra bonds from text such as ``"0-2, 0-3"`` or ``"[[0, 2], [0, 3]]"``.

    Pairs of equal indices are dropped, since an atom is never bonded to itself.
    """
    pairs: list[tuple[int, int]] = []
    for first, second in BOND_PAIR_PATTERN.findall(str(text)):
        pair = (int(first), int(second))
        if pair[0] != pair[1] and pair not in pairs:
            pairs.append(pair)
    return tuple(pairs)


def _parse_flag(text: str, default: bool) -> bool:
    value = text.strip().lower()
    if value in _TRUE:
        return True
    if value in _FALSE:
        return False
    return default


def _parse_positive_integer(name: str, text: str, default: int) -> int:
    try:
        value = int(text.strip())
    except (AttributeError, ValueError):
        return default
    if value < 1:
        raise ValueError(f"{name} must be a positive integer, got {text!r}")
    return value


@dataclass(frozen=True)
class LoaderSettings:
    """Everything that defines one load.

    Args:
        filepath: Structure file to read; any format supported by ASE works.
        explicit_bonds: Extra bonds to draw, as ``(atom, atom)`` pairs of Python
            indices (the first atom of the file is 0). Bond detection is based on
            covalent radii, so it misses some bonds of metals and of high-energy
            geometries.
        slow_factor: Insert a keyframe every ``slow_factor`` frames.
        include_first_snapshot: Keyframe every snapshot of the file. When
            ``False`` the first snapshot is only used as the static geometry of
            the scene, which is how the loader behaved before.
        collection_name: Blender collection that holds the new objects; the name
            of the file is used when empty.
    """

    filepath: Path
    explicit_bonds: tuple[tuple[int, int], ...] = ()
    slow_factor: int = 1
    include_first_snapshot: bool = True
    collection_name: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "filepath", Path(self.filepath).expanduser())
        object.__setattr__(
            self,
            "explicit_bonds",
            tuple((int(first), int(second)) for first, second in self.explicit_bonds),
        )
        if self.slow_factor < 1:
            raise ValueError(f"slow_factor must be a positive integer, got {self.slow_factor}")

    @property
    def collection(self) -> str:
        """Name of the collection that will be created for the structure."""
        return self.collection_name.strip() or self.filepath.stem or "Molecule"

    @classmethod
    def from_environment(cls, filepath=None, **overrides) -> LoaderSettings:
        """Build the settings from the ``MOLECULAR_LOADER_*`` environment variables.

        Arguments win over the environment, which wins over the defaults.
        """
        environment = os.environ
        values = {
            "filepath": filepath or environment.get(FILE_VARIABLE, ""),
            "explicit_bonds": parse_bond_pairs(environment.get(EXPLICIT_BONDS_VARIABLE, "")),
            "slow_factor": _parse_positive_integer(
                SLOW_FACTOR_VARIABLE, environment.get(SLOW_FACTOR_VARIABLE, ""), 1
            ),
            "include_first_snapshot": _parse_flag(
                environment.get(INCLUDE_FIRST_VARIABLE, ""), default=True
            ),
            "collection_name": environment.get(COLLECTION_VARIABLE, ""),
        }
        values.update(overrides)
        if not values["filepath"]:
            raise ValueError(
                "no structure file given: pass a path or set the "
                f"{FILE_VARIABLE} environment variable"
            )
        return cls(**values)
