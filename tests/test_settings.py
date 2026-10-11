"""Tests of the option parsing, without Blender and without ASE."""

from __future__ import annotations

from pathlib import Path

import pytest

from molecular_loader.settings import LoaderSettings, parse_bond_pairs


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("", ()),
        ("0-2, 0-3", ((0, 2), (0, 3))),
        ("[[0, 2], [0, 3]]", ((0, 2), (0, 3))),
        ("0,2;0,3", ((0, 2), (0, 3))),
        ("0 2 3 4", ((0, 2), (3, 4))),
        ("7-9", ((7, 9),)),
    ],
)
def test_parse_bond_pairs(text, expected):
    assert parse_bond_pairs(text) == expected


def test_parse_bond_pairs_drops_self_bonds():
    assert parse_bond_pairs("2-2, 1-2") == ((1, 2),)


def test_collection_defaults_to_the_file_name():
    settings = LoaderSettings(filepath="some/path/ethanol.xyz")
    assert settings.collection == "ethanol"


def test_collection_can_be_given():
    settings = LoaderSettings(filepath="ethanol.xyz", collection_name="  solute  ")
    assert settings.collection == "solute"


def test_file_path_is_expanded_and_kept_as_path():
    settings = LoaderSettings(filepath="~/ethanol.xyz")
    assert settings.filepath == Path.home() / "ethanol.xyz"


def test_slow_factor_must_be_positive():
    with pytest.raises(ValueError, match="slow_factor"):
        LoaderSettings(filepath="ethanol.xyz", slow_factor=0)


def test_explicit_bonds_are_normalised():
    settings = LoaderSettings(filepath="ethanol.xyz", explicit_bonds=[["1", "2"]])
    assert settings.explicit_bonds == ((1, 2),)


def test_from_environment(monkeypatch):
    monkeypatch.setenv("MOLECULAR_LOADER_FILE", "/tmp/molecule.xyz")
    monkeypatch.setenv("MOLECULAR_LOADER_SLOW_FACTOR", "4")
    monkeypatch.setenv("MOLECULAR_LOADER_EXPLICIT_BONDS", "0-2, 0-3")
    monkeypatch.setenv("MOLECULAR_LOADER_COLLECTION", "solute")
    monkeypatch.setenv("MOLECULAR_LOADER_INCLUDE_FIRST_SNAPSHOT", "no")

    settings = LoaderSettings.from_environment()

    assert settings.filepath == Path("/tmp/molecule.xyz")
    assert settings.slow_factor == 4
    assert settings.explicit_bonds == ((0, 2), (0, 3))
    assert settings.collection == "solute"
    assert settings.include_first_snapshot is False


def test_from_environment_arguments_win(monkeypatch):
    monkeypatch.setenv("MOLECULAR_LOADER_FILE", "/tmp/from-env.xyz")
    monkeypatch.setenv("MOLECULAR_LOADER_SLOW_FACTOR", "4")

    settings = LoaderSettings.from_environment(filepath="/tmp/from-argument.xyz", slow_factor=1)

    assert settings.filepath == Path("/tmp/from-argument.xyz")
    assert settings.slow_factor == 1


def test_from_environment_needs_a_file(monkeypatch):
    monkeypatch.delenv("MOLECULAR_LOADER_FILE", raising=False)
    with pytest.raises(ValueError, match="no structure file"):
        LoaderSettings.from_environment()


def test_from_environment_ignores_garbage(monkeypatch):
    monkeypatch.setenv("MOLECULAR_LOADER_FILE", "/tmp/molecule.xyz")
    monkeypatch.setenv("MOLECULAR_LOADER_SLOW_FACTOR", "many")
    monkeypatch.setenv("MOLECULAR_LOADER_INCLUDE_FIRST_SNAPSHOT", "maybe")

    settings = LoaderSettings.from_environment()

    assert settings.slow_factor == 1
    assert settings.include_first_snapshot is True
