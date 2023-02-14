# Installation 

```
pip install callingcardstools
```

To start using the command line tools, see the help message with:

```
callingcardstools --help
```

# Development Installation

1. install [poetry](https://python-poetry.org/)
  - I prefer to set the default location of the virtual environment to the 
  project directory. You can set that as a global configuration for your 
  poetry installation like so: `poetry config virtualenvs.in-project true`

2. git clone the repo

3. cd into the repo and issue the command `poetry install`

4. shell into the virtual environment with `poetry shell`

5. build the package with `poetry build`

6. install the callingcardstools packge into your virtual environment 
  `pip install dist/callingcardstools-...`
  - Note: you could figure out how to use the pip install `-e` flag to 
  have an interactive development environment. I don't think that is compatible 
  with only the `pyproject.toml` file, but if you look it up, you'll find good 
  stackoverflow instructions on how to put a dummy `setup.py` file in to make 
  this possible

7. Building the Dockerimage:

Currently the Dockerimage is built from a stable version on github

Note that unless I set it up, you won't be able to push to my dockerhub repo. 
I think that is possible to do, though. If you wish to push to your own dockerhub, 
replace the cmatkhan to your username.

```bash
docker build -t cmatkhan/callingcardstools - < Dockerfile
```

where cmatkhan/callingcardstools is the tag. This will default to the version 
`latest`

To push:

```bash
docker push cmatkhan/callingcardstools
```