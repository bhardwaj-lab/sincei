
## sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

[![DOI](https://zenodo.org/badge/271841139.svg)](https://zenodo.org/badge/latestdoi/271841139) [![Documentation Status](https://readthedocs.org/projects/sincei/badge/?version=latest)](https://sincei.readthedocs.io/en/latest/?badge=latest) [![PyPI Version](https://img.shields.io/pypi/v/sincei.svg?style=plastic)](https://pypi.org/project/sincei/) [![test](https://github.com/vivekbhr/sincei/actions/workflows/test.yml/badge.svg)](https://github.com/vivekbhr/sincei/actions/workflows/test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## [Full Documentation](http://sincei.rtfd.io/).

## Installation

sincei is a command line toolkit based on python3, and can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Create a new conda environment and install sincei stable release from github using:

```
conda create -n sincei -c anaconda python=3.8
conda activate sincei
(sincei): pip install --editable=git+https://github.com/vivekbhr/sincei.git@master#egg=sincei
```

For the development version, try:

```
(sincei): pip install --editable=git+https://github.com/vivekbhr/sincei.git@develop#egg=sincei
```

## Usage

**Get the tool list with `sincei --help`**

Each tool begins with the prefix sc<tool_name>, such as:

 $ scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage


