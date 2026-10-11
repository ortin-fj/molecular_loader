# molecular_loader

Welcome to Molecular Loader for Blender! Blender already has a built-in function
to upload molecules, but this one can also upload molecular dynamics
trajectories, so you can make beautiful movies for your scientific conferences.
You can load any type of structure file that ASE supports (XYZ, PDB, CIF...).

https://github.com/user-attachments/assets/e9e372c9-1ff0-4b95-8ca8-a1e5bcb6e2fc

The package "Molecular Nodes" from Brady Johnston
(https://github.com/BradyAJohnston/MolecularNodes) is far superior in this
aspect, but this code allows to easily have full control of every single bond in
the simulation.

I would like to acknowledge "Blender Atomic Loader"
(https://github.com/nanotech-empa/blender-atomic-loader). My code is more of an
update of his code than a complete project of mine. His code was designed to load
single frames in Blender v.2 and unfortunately it doesn't work anymore (at least
in my PC). My code updates his project so it works in newer versions of Blender,
but it also adds the functionality to read molecular dynamics trajectories and
other minor changes.


## Installation

### 1. Get the add-on

Download the zip of a release, or build it from the sources:

```bash
uv run python scripts/build_addon.py      # writes dist/molecular_loader-<version>.zip
```

Then, in Blender: **Edit > Preferences > Add-ons > Install from Disk** and pick
the zip. The same archive works as an extension (Blender 4.2 and newer) and as a
legacy add-on (Blender 3.6 - 4.1); `scripts/build_addon.py` validates it with
`blender --command extension validate`.

### 2. Give Blender the ASE library

Blender brings its own Python interpreter, which is *not* the interpreter of this
project, so ASE has to be installed inside it. Anything supported by ASE (XYZ,
PDB...) needs both `ase` and the `scipy` that ASE imports.

The easy way is the **Install dependencies** button of the Molecular Loader panel
(3D viewport sidebar, <kbd>N</kbd>), which installs them into the `libs`
directory of the add-on. From the sources it is the same thing:

```bash
uv run python scripts/install_dependencies.py
```

Or manually, with any Blender:

```bash
uv pip install --python "$(blender --background --python-expr 'import sys; print(sys.executable)')" \
    --target <add-on directory>/libs --no-deps ase scipy
```

The `--no-deps` is deliberate: `ase` and `scipy` are installed without a copy of
numpy, so they use the numpy that Blender ships. Without a matching numpy inside
Blender's interpreter, the compiled parts of scipy would not import.


## Usage

### From Blender: `File > Import > Molecular Structure`

The file browser opens, and the structure is built in the scene. The options of
the panel are used as defaults, and the panel remembers the last file.

### From the panel

Open the sidebar of the 3D viewport (<kbd>N</kbd>) and look for the *Molecular
Loader* tab: it has the file path, the options below, a **Choose file...** button
(browser), a **Load** button (uses the path of the panel) and the state of the
libraries. Nothing has to be edited by hand.

### From the command line

```bash
blender --background --python scripts/load_in_blender.py -- molecule.xyz
blender --background --python scripts/load_in_blender.py -- molecule.xyz \
        --slow-factor 5 --explicit-bonds "0-2, 0-3" --collection solute \
        --save-blend scene.blend
```

Every option can also come from the environment, which makes it easy to drive the
loader from another program:

```bash
MOLECULAR_LOADER_FILE=molecule.xyz MOLECULAR_LOADER_SLOW_FACTOR=5 \
    blender --background --python scripts/load_in_blender.py
```

### From a script (or Blender's Text Editor)

With the add-on installed, its package is importable from any Python script that
runs inside Blender:

```python
import molecular_loader

molecular_loader.load("molecule.xyz")                            # defaults
molecular_loader.load("molecule.xyz", slow_factor=5, collection_name="solute")
report = molecular_loader.load_from_environment()                # MOLECULAR_LOADER_* variables
print(report.summary())
```


## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `filepath` | - | Structure file; any format supported by ASE (xyz, pdb, cif...) |
| `explicit_bonds` | `()` | Bonds that are not detected, as pairs of atom indices: `[(0, 2), (0, 3)]` in Python, `"0-2, 0-3"` in the panel or the command line |
| `slow_factor` | `1` | Insert a keyframe every N frames |
| `include_first_snapshot` | `True` | Keyframe the first snapshot as well. With `False` the first snapshot is only the static geometry, which is what the old script did |
| `collection_name` | file name | Collection of the new objects |

The same names are available as environment variables:
`MOLECULAR_LOADER_FILE`, `MOLECULAR_LOADER_EXPLICIT_BONDS`,
`MOLECULAR_LOADER_SLOW_FACTOR`, `MOLECULAR_LOADER_INCLUDE_FIRST_SNAPSHOT`,
`MOLECULAR_LOADER_COLLECTION` and `MOLECULAR_LOADER_LIBS` (extra directories with
`ase` and `scipy`).


## What ends up in the scene

* The objects go into a new collection, named after the file (or after
  `collection_name`). A new collection is always created, so an existing scene is
  never modified.
* Object names follow the atoms they represent: `Atom1`, `Atom2`... for the
  spheres and `bond1_2` for the half of the bond between atoms 1 and 2 that is
  coloured like atom 1. A bond between two different elements is drawn as two
  half-cylinders (`bond1_2` and `bond2_1`), while a bond between two equal
  elements is a single cylinder.
* Bonds are detected from the covalent radii of ASE. For high-energy structures
  or species containing metals not all the bonds may be drawn; list the missing
  ones in `explicit_bonds`. Indices are Python indices, so the first atom of the
  file is 0.
* Each element gets one Principled BSDF material (`C_mat`, `H_mat`...). Elements
  that are not in the table of the package get a neutral grey fallback.
* Objects of the same element share their mesh, and so do the cylinders of the
  same element and radius, so even big trajectories only need a handful of
  meshes.
* One keyframe per snapshot is inserted, so one snapshot is shown per frame of
  the movie (or every `slow_factor` frames). The frame range of the scene grows
  to cover the animation, but it is never shrunk.
* A file with a single geometry gets no keyframes at all.
* The scene is switched to Cycles. This code is based on a brute-force approach:
  large files (I got to upload 700 atoms in 500 snapshots, i.e. thousands of
  objects in Blender) take a while to load.
* Once the keyframes have been assigned, it is annoying to change the orientation
  of the species. Make sure to rotate the species in your structure file before
  creating the final scene.


## Development

The project is managed with [uv](https://docs.astral.sh/uv/):

```bash
uv sync --extra dev          # virtual environment with ASE, pytest and ruff
uv run pytest                # options, bond detection and geometry (no Blender)
uv run ruff check .          # lint
uv run python tests/data/make_fixtures.py           # regenerate the test files
```

The tests that need Blender run *inside* Blender and build a real scene:

```bash
MOLECULAR_LOADER_LIBS=src/molecular_loader/libs \
    blender --background --factory-startup --python tests/blender/smoke_test.py
```

### Layout

```
blender_manifest.toml              Metadata of the add-on (Blender 4.2+ extensions)
pyproject.toml                     Project metadata, dependencies and tool settings
scripts/build_addon.py             Build (and validate) the installable zip
scripts/install_dependencies.py    Install ASE where Blender can import it
scripts/load_in_blender.py         Load a structure headless, from the command line
src/molecular_loader/
    __init__.py                    Public API, and the register/unregister of the add-on
    ui.py                          File > Import, the sidebar panel and its operators
    loader.py                      Builds the whole scene of a structure (entry point)
    settings.py                    LoaderSettings: files, extra bonds, environment
    structure.py                   Reads structures and trajectories with ASE
    bonds.py                       Decides which pairs of atoms are bonded
    geometry.py                    Where a bond cylinder goes and how it is rotated
    builder.py                     Spheres, cylinders, shared meshes and the collection
    animation.py                   The keyframes of a trajectory
    materials.py                   One Principled BSDF material per element
    elements.py                    Colour and radius of every element
    deps.py                        Finds (and installs) ASE and scipy inside Blender
tests/                             pytest tests, and the Blender smoke test in tests/blender
```

### uv and Blender: what is compatible and what is not

`bpy`, `bmesh` and `mathutils` only exist inside Blender, so **the virtual
environment is for tools and tests, never for running the loader**: `uv run
python -c "import bpy"` cannot work, and the loader is always executed by
Blender. The version of the interpreter of the virtual environment (3.12 here)
does not have to match the one of Blender (3.13 in Blender 5.2, 3.10 in Blender
3.6) either, because the loader never uses the virtual environment.

What `uv` is good for here:

* `uv sync`, `uv run pytest` and `uv run ruff`: development, with a Python that
  is independent of the installed Blender.
* `uv pip install --python <blender python> --target <libs>`: the recommended way
  to install ASE *for Blender*. `uv` reads the version and the platform of the
  interpreter given by `--python`, so the wheels match that Blender, and it
  reuses them from its cache instead of downloading them again.

`scripts/install_dependencies.py` does that (`--blender` to point at a specific
installation, `--python` to skip the search, `--no-uv` to fall back to
`python -m pip`).


## Troubleshooting

* **"ASE is not available inside Blender's Python"**: install the dependencies as
  explained above; the button in the panel does it from inside Blender.
* **The panel says "Missing libraries: ase, scipy"**: the same thing, press
  *Install dependencies* (it needs internet).
* **Nothing happens when running a script**: one snapshot is keyframed per frame,
  so a file with a single geometry has no animation at all.
* **Blender cannot import scipy**: the installed scipy needs a newer numpy than
  the one Blender ships. Update Blender, or install an older scipy with
  `uv pip install --python <blender python> --target <libs> --no-deps "scipy<1.16"`.
