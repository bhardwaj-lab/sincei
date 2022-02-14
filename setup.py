# -*- coding: utf-8 -*-

import re

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


def openREADME():
    """
    This is only needed because README.md is UTF-8 encoded and that won't work
    under python3 iff sys.getfilesystemencoding() returns 'ascii'
    Since open() doesn't accept an encoding in python2...
    """
    try:
        f = open("README.md", encoding="utf-8")
    except:
        f = open("README.md")

    foo = f.read()
    f.close()
    return foo


setup(
    name='sincei',
    version=get_version(),
    author=' Vivek Bhardwaj',
    author_email='vivbhr@gmail.com',
    packages=find_packages(),
    scripts=['bin/scBulkCoverage.py', 'bin/scClusterCells.py', 'bin/scCountQC.py',
             'bin/scCountReads.py', 'bin/scFilterBarcodes.py', 'bin/scFilterStats.py',
             ],
    include_package_data=False,
    url='https://github.com/vivekbhr/sincei',
    license='LICENSE.txt',
    description='A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data ',
    long_description=openREADME(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        "scanpy",
        "deeptools",
        "gensim"
    ],
    zip_safe=True,
    cmdclass={'sdist': sdist, 'install': install}
)
