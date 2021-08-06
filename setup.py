from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install


# Classes for shortcut executable/application
class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        develop.run(self)
        from pyshortcuts import make_shortcut
        make_shortcut(
            script='gui/evSeq_gui.py',
            name='evSeq',
        )

class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        install.run(self)
        from pyshortcuts import make_shortcut
        make_shortcut(
            script='gui/evSeq_gui.py',
            name='evSeq',
        )

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

with open('evSeq/__init__.py', 'r') as f:
    init = f.readlines()

for line in init:
    if '__author__' in line:
        __author__ = line.split("'")[-2]
    if '__email__' in line:
        __email__ = line.split("'")[-2]
    if '__version__' in line:
        __version__ = line.split("'")[-2]

setup(
    name='evSeq',
    version=__version__,
    author=__author__,
    author_email=__email__,
    description='Package for analyzing pooled NGS sequencing results.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'colorcet',
        'holoviews',
        'bokeh',
        'numpy',
        'pandas',
        'tqdm',
        'scipy',
        'gooey',
        'jupyterlab',
        'pyshortcuts',
        'ninetysix'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'evSeq = evSeq.cmd:main'
        ],
        'gui_scripts': [
            'evSeq-GUI = gui.evSeq_gui:main'
        ]
    },
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    }
)
