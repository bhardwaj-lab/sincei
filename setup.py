# -*- coding: utf-8 -*-

import re
# read the contents of your README file
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.
__version__ = '%s'
"""


def get_version():
    try:
        f = open("sincei/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)


class install(_install):

    def run(self):
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return


this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='sincei',
    version=get_version(),
    author=' Vivek Bhardwaj',
    author_email='vivbhr@gmail.com',
    packages=find_packages(),
    scripts=['bin/scBulkCoverage', 'bin/scClusterCells', 'bin/scCountQC',
             'bin/scCountReads', 'bin/scFilterBarcodes', 'bin/scFilterStats',
             ],
    include_package_data=False,
    url='https://github.com/vivekbhr/sincei',
    license='LICENSE.txt',
    description='A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data ',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        "scanpy",
        "loompy",
        "deeptools",
        "gensim",
        "leidenalg",
        "networkx",
        "community",
        "igraph"
    ],
    zip_safe=True,
    cmdclass={'sdist': sdist, 'install': install}
)
