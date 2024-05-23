
<img align="right" src="./docs/content/images/sincei-logo.png">


## sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

[![DOI](https://zenodo.org/badge/271841139.svg)](https://zenodo.org/badge/latestdoi/271841139) [![Documentation Status](https://readthedocs.org/projects/sincei/badge/?version=latest)](https://sincei.readthedocs.io/en/latest/?badge=latest) [![PyPI Version](https://img.shields.io/pypi/v/sincei.svg?style=plastic)](https://pypi.org/project/sincei/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![test](https://github.com/vivekbhr/sincei/actions/workflows/test.yml/badge.svg)](https://github.com/vivekbhr/sincei/actions/workflows/test.yml) [![codecov](https://codecov.io/gh/vivekbhr/sincei/graph/badge.svg?token=VRTMITHHBI)](https://codecov.io/gh/vivekbhr/sincei)

## Features

sincei provides a flexible, easy-to-use command-line interface to work with single-cell data directly from BAM files. It can:

 - Aggregate signal in bins, genes or any feature of interest from single-cells.
 - Perform read-level and count-level quality control.
 - Perform dimentionality reduction and clustering of all kinds of single-cell data (open chromatin, histone marks, methylation, gene expression etc..).
 - Create coverage files (bigwigs) for visualization.

## [Full Documentation](http://sincei.rtfd.io/)

Please browse the full documentation for tutorials on how to use sincei on command line, as well as details of our python API.

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

```
scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage
```

## Questions and discussions

To ask  a question related to sincei or start a new discussion, please use our [github discussion forum](https://github.com/vivekbhr/sincei/discussions).
