# Installation
For non-programmers and those unfamiliar with [Anaconda](https://www.anaconda.com/) or [GitHub](https://www.github.com), see the [programing basics page](basics.md) for information on how to set up your computer environment to run `evSeq`.
## Installing from GitHub with the `conda` environment
The recommended way to install `evSeq` invloves cloning the reopsitory from GitHub and then creating the [`conda` environment](https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) from the included `evSeq.yml` file, as follows:
```
git clone https://github.com/fhalab/evSeq.git
cd evSeq
conda env create -f evSeq.yml
```

### Using the `evSeq` environment
This environment can subsequently be activated any time you want to run `evSeq` from the command line with:
```
conda activate evSeq
```
and deactivated with
```
conda deactivate
```
## Standard `pip` Install and Dependencies
Advanced users: If you would rather not use the `evSeq` environment described above and run in a custom environment (or, if you're a brave soul, your base environment), below are the `evSeq` dependencies, most of which are available through `conda` (only `ninetysix` and `pyshortcuts` require `pip`, but they are small packages with few dependencies).

The `evSeq` dependencies are explicitly listed in the `evSeq.yml` environment file:
```yml
# evSeq without GUI
- biopython>=1.78
- colorcet
- holoviews>=1.14.3
- bokeh>=2.3.2
- numpy
- pandas
- python>=3.7
- tqdm
- scipy
- ninetysix

# With GUI, also need:
- gooey
- pyshortcuts
```

Releases of `evSeq` are hosted on the [Python Package Index (PyPI)](https://pypi.org/project/evseq/), and thus can be `pip` installed in the standard way. From any location, run
```
pip install evSeq
```

## Confirming your installation
Once `evSeq` is installed and your environment is suitable, you should [confirm your installation](comp/usage.md#confirming-your-installation).

---

*Next page: [Running `evSeq`](usage.md).*

*Back to the [main page](../index.md).*