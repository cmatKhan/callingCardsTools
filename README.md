# callingCardsTools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10042067.svg)](https://doi.org/10.5281/zenodo.10042067)
[![callingCardsTools Test Coverage](https://github.com/cmatKhan/callingCardsTools/actions/workflows/codecov.yml/badge.svg?branch=main)](https://github.com/cmatKhan/callingCardsTools/actions/workflows/codecov.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/callingcardstools/README.html)

## Introduction

`CallingCardsTools` Provides both an API and a number of cmd line tools
for processing raw Calling Cards data. This is used in the
[nf-core/callingcards](https://github.com/nf-core/callingcards) pipeline,
which provides a workflow to process both yeast and mammals Calling Cards data.

## Documentation

[Served Documentation](https://cmatkhan.github.io/callingCardsTools/) provides
information on filetypes and the API. For help with the cmd line tools,
simply install callingcardstools (see below) and do:

```bash
callingcardstools --help
```

Each of the cmd line tools also provides a `--help` message.

## Installation

callingCardsTools is available through bioconda:

```bash
conda install -c bioconda callingcardstools
```

pypi:

```bash
pip install callingcardstools
```

or github (this will be the most current version):

```bash
pip install git+https://github.com/cmatkhan/callingCardsTools.git
```

After installing, you can get help with the cmd line tools by doing:

```bash
callingcardstools --help
```

## Callingcardstools is containerized:

- A singularity container is hosted on
  [Galaxyhub](https://depot.galaxyproject.org/singularity/). If you go to this
  site, make sure the 'c's have loaded and then search for 'callingcardstools'.
  There is a container for each version which is on bioconda. Make sure you get
  the correct version.

- A docker container is hosted on
  [quay (and biocontainers)](https://quay.io/repository/biocontainers/callingcardstools).
  Again, make sure you get the correct version.

## Development Installation

1. install [poetry](https://python-poetry.org/)

- I prefer to set the default location of the virtual environment to the
  project directory. You can set that as a global configuration for your
  poetry installation like so: `poetry config virtualenvs.in-project true`

2. git clone the repo

3. cd into the repo and issue the command `poetry install`

4. shell into the virtual environment with `poetry shell`

5. you can `pip install -e .` to install the package in editable mode. This is
   useful if you want to test the cmd line interface as you make changes to the
   source code.
