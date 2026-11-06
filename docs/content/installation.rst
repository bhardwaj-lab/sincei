Installation
============

Sincei is a command line toolkit written in Python and Rust. The stable version of
sincei is distributed as a Python binary wheel and can be installed using
`conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ ,
`pip <https://packaging.python.org/en/latest/tutorials/installing-packages/#installing-from-pypi>`_,
or `uv <https://docs.astral.sh/uv/>`_.
Development versions can be installed and compiled from github using
`maturin <https://github.com/pyo3/maturin>`_.

conda
-----

Bioconda hosts the latest stable version of sincei. Create a new conda environment
and install sincei in it using:

.. code-block:: bash

    conda create -n sincei -c bioconda -c conda-forge sincei

PyPi
----

You can also install the latest stable version of sincei from PyPi using pip or uv:

.. code-block:: bash

    pip install sincei
    # or
    uv pip install sincei

Development version
-------------------

You can install the development version of sincei by cloning our
`github repository <https://github.com/bhardwaj-lab/sincei>`_ and compiling it using
maturin. This will install the latest development version of sincei in your conda
environment:

.. code-block:: bash

    conda create -n sincei python=3.12 uv maturin
    conda activate sincei
    git clone -b rust_release https://github.com/bhardwaj-lab/sincei.git
    cd sincei
    maturin develop --release --uv
