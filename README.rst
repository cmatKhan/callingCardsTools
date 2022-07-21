# System Dependencies

[Bedtools](https://github.com/arq5x/bedtools2). 
Download and install following their instructions. It must be 
available in your `$PATH`

# General use

Once you have bedtools in your path, you can install like so:

```
pip install callingcardstools
```

There is currently one entry point (command line script):

```
macs_peak_caller
```

This will be available in your `$PATH` after you install the package.

Example usage:

```
macs_peak_caller \
    -e human/TAG_AY53-1_50k_downsampled_human_map_sort_SUBSET.ccf \
    -t ttaa_subset_99_rows.tsv \
    -a ref_subset_500_plus_chrm.bed \
    -b TAG_background_SUBSET.ccf
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