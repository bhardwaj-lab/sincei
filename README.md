
## sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

[![DOI](https://zenodo.org/badge/271841139.svg)](https://zenodo.org/badge/latestdoi/271841139) [![Documentation Status](https://readthedocs.org/projects/sincei/badge/?version=latest)](https://sincei.readthedocs.io/en/latest/?badge=latest) [![test](https://github.com/vivekbhr/sincei/actions/workflows/test.yml/badge.svg)](https://github.com/vivekbhr/sincei/actions/workflows/test.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## [Full Documentation](http://sincei.rtfd.io/).

## Installation

sincei is a command line toolkit based on python3, and can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Create a new conda environment and install sincei stable release from github using:

```
cd <your_programs_folder>
conda create -n sincei -c conda-forge -c bioconda scanpy deeptools
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


