<p align="center">
  <img
    src="https://raw.githubusercontent.com/moltools/RetroMol/main/logo.png"
    height="150"
    alt="RetroMol logo"
  >
</p>

<h1 align="center">
  RetroMol
</h1>

<p align="center">
    <a href="https://github.com/MolTools/RetroMol/actions/workflows/tests.yml">
      <img alt="testing & quality" src="https://github.com/MolTools/RetroMol/actions/workflows/tests.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/retromol">
      <img alt="PyPI" src="https://img.shields.io/pypi/v/retromol" /></a>
    <a href="https://pypi.org/project/retromol">
      <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/retromol" /></a>
     <a href="https://doi.org/10.5281/zenodo.17410940">
      <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17410940.svg" alt="DOI" /></a>
</p>

RetroMol is retrosynthetic parsing and fingerprinting tool for modular natural products. 

RetroMol is designed to facilitate clustering modular natural products based on biosynthetic similarity, and to enable cross-modal retrieval between modular natural products and their coding biosynthetic gene clusters.

RetroMol powers the RetroMol webapp, available at [here](https://retromol.bioinformatics.nl/). The webapp can be used to query NPAtlas, MIBiG, and antiSMASH database entries for modular natural products based on biosynthetic similarity.

## Installation

The most recent code and data can be installed directly from GitHub with:

```shell
pip install git+https://github.com/MolTools/RetroMol.git
```

The latest stable release can be installed from PyPI with:

```shell
pip install retromol
```

RetroMol has been developed for Linux and MacOS.

## Getting started

The retromol command line tool is automatically installed alongside installation of the package. The command line tool can be used from the shell with the `--help` flag to show all subcommands:

```shell
python3 -m retromol --help
```

The RetroMol CLI has two modes, single and batch:
* `single`: process a single compound at a time.
* `batch`: process multiple compounds in a single command.

In either case a the output folder will contain a log file together with the results in either JSON or JSONL format. JSONL is standard output mode for batch mode to allow for easy parsing of large result sets. Batch mode also supports parallel processing.

Any column, field, or property in the input file, either CSV, TSV, SDF, or JSON, is preserved as props in the output JSON or JSONL.

Stereochemistry parsing is supported by supplying the `--matchstereochem` flag. This flag annotates the identifier of every identified monomer with R/S and E/Z annotation where applicable.

Result JSONs or lines from a JSONL file can be loaded into Python using RetroMol's `Result` class for further downstream analyses:

```python
from retromol.io import Result

result = Result.from_serialized(<json_dict>)

# e.g., calculate coverage
coverage = result.best_total_coverage()
```

Check out the [examples](https://github.com/moltools/RetroMol/tree/main/examples) folder for example scripts demonstrating how to use RetroMol as a library:
* [Read out and align linear monomer readouts](https://github.com/moltools/RetroMol/tree/main/examples/align_compounds.py)
* [Calculate and cluster monomer fingerprints](https://github.com/moltools/RetroMol/tree/main/examples/cluster_compounds.py)

### Using custom rules

RetroMol comes with a default set of retrosynthetic rules for modular natural products.

You can also provide your own custom rules and configurations. See [the custom rules documentation](docs/customize_rules.md) for details on the YAML formats supported.

## Attribution

### License

The code in this package is licensed under the MIT License.

## For Developers

The final section of the README is for if you want to get involved by making a code contribution.

### Development installation

First fork the repository on GitHub, then clone your fork locally and install the package
in "editable" mode with the development dependencies:

```bash
git clone git+https://github.com/MolTools/RetroMol.git
cd RetroMol
pip install -e .[dev]
pip install hatch     # if you don't have hatch installed yet; needed for building the package
hatch env create dev  # create the development environment
```

You can now make code changes locally and have them immediately available for testing.

After testing your changes, you can commit and push them to your fork, and then open a pull request
on the main repository explaining your changes.

### Testing

After cloning the repository, the unit tests in the `tests/` folder can be run
reproducibly with:

```shell
pytest tests
```

Additionally, these tests are automatically re-run with each push and pull request on `main` in a
[GitHub Action](https://github.com/MolTools/RetroMol/actions?query=workflow%3ATests).
