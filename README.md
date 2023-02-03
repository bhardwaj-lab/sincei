
<img src="./docs/images/sincei.png" alt="sincei logo" style="height: 200px; width:400px;"/>

## Bhardwaj V. (2022) sincei: A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data.

[![DOI](https://zenodo.org/badge/271841139.svg)](https://zenodo.org/badge/latestdoi/271841139)

 - **Documentation**: [![Documentation Status](https://readthedocs.org/projects/sincei/badge/?version=latest)](https://sincei.readthedocs.io/en/latest/?badge=latest)
 - **Source code**: [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Installation

sincei is a command line toolkit based on python3, and can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Create a new conda environment and install sincei using:

```
cd <programs_folder>
conda create -n sincei -c bioconda -c conda-forge scanpy deeptools
conda activate sincei
(sincei): pip install --editable=git+https://github.com/vivekbhr/sincei.git@master#egg=sincei
```

## Usage

**Get the tool list with `sincei --help`**

Each tool begins with the prefix sc<tool_name>, such as:

 $ scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage

Read the [full documentation here]().
