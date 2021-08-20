# Computational Basics
## Useful information for non-programmers

This section details installation of high level dependencies: `Git Bash` (Windows users), `git`, and `Anaconda`. If you have installed and are familiar with these items, you can skip this section and move on to [Installation](3-installation.md). Installation on Linux is not detailed here, as we just assume you know what you're doing. Once this section is completed, continue to [installation](3-installation.md).

### Installing `Git`
#### Windows users
Install `Git`/`Git Bash` by following the instructions [here](https://www.stanleyulili.com/git/how-to-install-git-bash-on-windows/). Unless you know what you're doing, sticking to the default options during install is best.

#### Mac users
Install `Git` by following the instructions [here](https://www.atlassian.com/git/tutorials/install-git), if not already installed on your machine. Some commands must be performed through the terminal. You can open a terminal by typing command + spacebar to open Spotlight Search, then type "Terminal"; opening the search result will give you a terminal window.

### Opening a Terminal Window
Later steps in our installation will require you to write commands in terminal. To open a terminal in Windows, right-click on either the desktop background or within your file explorer, then open `Git Bash`. For Mac users, type command + spacebar, search "Terminal", then open the search result.

### Installing Anaconda
`Anaconda` is an open-source package managment framework for scientific computing with Python. For details, look at their website [here](https://www.anaconda.com/). All software that supports `evSeq` is or can be handled by the `Anaconda` package manager. See below for installation instructions on both Windows and Mac.

#### Windows users
Install `Anaconda` following the instructions [here](https://docs.anaconda.com/anaconda/install/windows/). At step 8, we recommend adding `Anaconda` to your `PATH` environment variable. Note that this is in contrast to the recommendation of `Anaconda`, but their concerns shouldn't apply for this use case. Once installation is complete, open a terminal window and enter the below command
```
conda init bash
```

#### Mac users
Install `Anaconda` following the instructions [here](https://docs.anaconda.com/anaconda/install/mac-os/).

### Construction of a Folder for holding `Git` repositories
The next step will be to install `evSeq`, as described in [Installation](3-installation.md). The recommended way involves cloning the `evSeq` repository from `GitHub`. If you have not worked with Git repos before, we recommend creating a folder where you can store all of them. Wherever seems reasonable to you (most likely your home directory), create a folder called `git` or `GitRepos` or whatever seems best for you. Navigate to this folder (`cd folder_name`) and perform installation from there.

---

*Next page: [Installation](3-installation.md).*

*Back to the [main page](index.md).*