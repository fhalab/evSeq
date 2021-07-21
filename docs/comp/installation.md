# Installation
For non-programmers and those unfamiliar with [Anaconda](https://www.anaconda.com/) or [GitHub](https://www.github.com), see the [programing basics page](basics.md) for information on how to set up your computer environment to run `evSeq`.
## Installing from GitHub
Clone and install the respository via
```
git clone https://github.com/fhalab/evSeq.git
cd evSeq
pip install .
```
Which will install the local version of the package cloned from GitHub.

## Using the `evSeq` environment
For reproducibility, the `evSeq` repository contains the environment file `evSeq.yml` which will create a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) with exact, OS-dependent versions of packages that `evSeq` requires to function properly.

If you [cloned the repository](#installing-from-github) and are in the `evSeq` directory, you may set up this environment with the following command:
```
conda env create -f evSeq.yml
```
Otherwise download the file [here](../../evSeq.yml) and follow the same command above, making sure that the path to the file's download location is included in the file argument (e.g., replace `evSeq.yml` with `Downloads/evSeq.yml` if the file is in the directory `Downloads`).

This environment can subsequently be activated any time you want to run `evSeq` with:
```
conda activate evSeq
```
and deactivated with
```
conda deactivate
```

## Installing from PyPI
#### (Standard `pip` install)

Releases of `evSeq` are also hosted on the [Python Package Index (PyPI)](https://pypi.org/project/evseq/), and thus can be `pip` installed in the standard way. From any location, run

```
pip install evSeq
```
#### Dependencies
Advanced users: If you would rather not use the `evSeq` environment described above and run in a custom environment (or, if you're a brave soul, your base environment), below are the `evSeq` dependencies, most of which are available through `conda` (only `ninetysix` requires `pip`). These dependencies are, of course, explicitly listed in the `evSeq.yml` environment file.

    # evSeq without GUI
    - biopython>=1.78
    - colorcet
    - holoviews
    - bokeh>=1.4.0
    - numpy
    - pandas
    - python>=3.7
    - tqdm
    - scipy
    - ninetysix

    # With GUI, also need:
    - gooey
## Confirming your installation
Once `evSeq` is installed and your environment is suitable, you should [confirm your installation](comp/usage.md#confirming-your-installation).

---

*Back to the [main page](../index.md).*