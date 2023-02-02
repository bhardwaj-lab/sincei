# -*- coding: utf-8 -*-

import re
# read the contents of your README file
#from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.
__version__ = '%s'
"""

VPATH="sincei/_version.py"
def get_version(path=VPATH):
    try:
        f = open(path)
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None

def get_dependencies():
    try:
        f = open("sincei/requirements.txt")
    except FileNotFoundError:
        return None
    out = []
    for line in f.readlines():
        out.append(line)
    return out

class sdist(_sdist):

    def run(self):
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)


class install(_install):

    def run(self):
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))

#this_directory = Path(__file__).parent
#long_description = (this_directory / "README.md").read_text()
long_description = open('README.md').read()

setup(
    name='sincei',
    version=get_version(),
    author=' Vivek Bhardwaj',
    author_email='vivbhr@gmail.com',
    packages=find_packages(),
    scripts=['bin/scBulkCoverage', 'bin/scClusterCells', 'bin/scCombineCounts', 'bin/scCountQC',
             'bin/scCountReads', 'bin/scFilterBarcodes', 'bin/scFilterStats', 'bin/scJSD', 'bin/sincei',
             ],
    include_package_data=True,
    package_dir={'sincei': 'sincei'},
    url='https://github.com/vivekbhr/sincei',
    license='LICENSE.txt',
    description='A user-friendly toolkit for QC, counting, clustering and plotting of single-cell (epi)genomics data ',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=get_dependencies(),
    setup_requires=get_dependencies(),
    zip_safe=False,
    cmdclass={'sdist': sdist, 'install': install}
)
