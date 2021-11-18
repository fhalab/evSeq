# Installation
For non-programmers and those unfamiliar with [Anaconda](https://www.anaconda.com/) or [GitHub](https://www.github.com), see the [programing basics page](2-basics.md) for information on how to set up your computer environment to run `evSeq`.
## Installing from GitHub with the conda environment
The recommended way to install `evSeq` invloves cloning the repository from GitHub and then creating the [`conda` environment](https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) from the included `evSeq.yml` file. Open a terminal window and navigate to the folder where you want to install `evSeq` (using the `cd` command to change directories as needed), then, in a terminal window, enter the following:
```
git clone https://github.com/fhalab/evSeq.git
cd evSeq
conda env create -f envs/evSeq.yml
```

This will create an environment with the correct dependencies and versions.

### Using the evSeq environment
This environment can subsequently be activated any time you want to run `evSeq` from the command line with:
```
conda activate evSeq
```
and deactivated with
```
conda deactivate
```
## Standard pip Install and Dependencies
Advanced users: If you would rather not use the `evSeq` environment described above and run in a custom environment (or, if you're a brave soul, your base environment), below are the `evSeq` dependencies, most of which are available through `conda` (only `ninetysix` and `pyshortcuts` require `pip`, but they are small packages with few dependencies).

The `evSeq` dependencies are explicitly listed all of the environment files and `setup.py` requirements:
```yml
- python>=3.7
- numpy
- pandas
- openpyxl
- biopython>=1.78
- scipy
- tqdm
- holoviews
- bokeh
- colorcet
- gooey
- pyshortcuts
- ninetysix
```
Jupyter installs are not listed here or in the `setup.py` requirements as they are not necessary if you do not plan to use any of the additional data visualization tools.

Releases of `evSeq` are hosted on the [Python Package Index (PyPI)](https://pypi.org/project/evseq/), and thus can be `pip` installed in the standard way. From any location, run
```
pip install evSeq
```

## Confirming your installation
Once `evSeq` is installed and your environment is suitable, you should [confirm your installation](4-usage.md#confirming-your-installation).

---

*Next page: [Running `evSeq`](4-usage.md).*

*Back to the [main page](index.md).*

# Updating
By default, `evSeq` is installed in non-dev mode. This means that changes to the code base on your computer will not be reflected come run-time. If you want an editable version of `evSeq`, install with the `evSeq_dev` environment (note, however, that this environment does not set exact versions of dependencies like `evSeq`). We recommend installing in non-dev mode (i.e., using the `evSeq` environment). To update `evSeq` when installed in a non-dev environment, the environment must be recreated. The below commands will update `evSeq`. First, navigate the evSeq repository folder via command line and enter the below commands:

```
git pull
conda remove -n evSeq --all
conda env create -f envs/evSeq.yml
```