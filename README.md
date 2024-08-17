
<img align="right" src="./docs/content/images/sincei-logo.png">


## sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

 [![Documentation Status](https://readthedocs.org/projects/sincei/badge/?version=latest)](https://sincei.readthedocs.io/en/latest/?badge=latest) [![PyPI Version](https://img.shields.io/pypi/v/sincei.svg?style=plastic)](https://pypi.org/project/sincei/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![test](https://github.com/vivekbhr/sincei/actions/workflows/test.yml/badge.svg)](https://github.com/vivekbhr/sincei/actions/workflows/test.yml) [![codecov](https://codecov.io/gh/vivekbhr/sincei/graph/badge.svg?token=VRTMITHHBI)](https://codecov.io/gh/vivekbhr/sincei) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sincei/README.html)

## Features

sincei provides a flexible, easy-to-use command-line interface to work with single-cell data directly from BAM files. It can:

 - Aggregate signal in bins, genes or any feature of interest from single-cells.
 - Perform read-level and count-level quality control.
 - Perform dimentionality reduction and clustering of all kinds of single-cell data (open chromatin, histone marks, methylation, gene expression etc..).
 - Create coverage files (bigwigs) for visualization.

For details, please [**read our preprint**](https://www.biorxiv.org/content/10.1101/2024.07.27.605424v1) describing sincei.

## [Full Documentation](http://sincei.rtfd.io/)

Please browse the full documentation for tutorials on how to use sincei on command line, as well as details of our python API.

## Installation

sincei is a command line toolkit based on python3, and can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

The recommended way to install sincei is via bioconda:

```
conda create -n sincei -c bioconda -c conda-forge sincei
```

Alternatively, a development version can be installed via GitHub.

```
conda create -n sincei -c anaconda python=3.8
conda activate sincei
(sincei): pip install --editable=git+https://github.com/bhardwaj-lab/sincei.git@master#egg=sincei
```

## Usage

**Get the tool list with `sincei --help`**

Each tool begins with the prefix sc<tool_name>, such as:

```
scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage
```

## Citation

Please cite sincei as: "Bhardwaj V. , Mourragui, S. (2024) User-friendly exploration of epigenomic data in single cells using sincei. biorXiv. doi: 10.1101/2024.07.27.605424"

## Questions and discussions

To ask  a question related to sincei or start a new discussion, please use our [github discussion forum](https://github.com/vivekbhr/sincei/discussions).
