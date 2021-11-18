from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install


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

class PostDevelopCommand(develop):
    """Post-installation for develop mode. Creates an executable/
    application to make a shortcut on the Desktop that opens the GUI
    and runs it from the install environment.
    """

    def run(self):
        develop.run(self)
        import build_shortcut
        build_shortcut.main()


class PostInstallCommand(install):
    """Post-installation for installation mode. Creates an executable/
    application to make a shortcut on the Desktop that opens the GUI
    and runs it from the install environment.
    """

    def run(self):
        install.run(self)
        import build_shortcut
        build_shortcut.main()

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
        'numpy',
        'pandas',
        'openpyxl',
        'biopython',
        'tqdm',
        'scipy',
        'holoviews',
        'bokeh',
        'colorcet',
        'gooey',
        'pyshortcuts',
        'ninetysix'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_data={'evSeq': ['util/index_map.csv']},
    entry_points={
        'console_scripts': [
            'evSeq = evSeq.cmd:main'
        ]
    },
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    }
)
