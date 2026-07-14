# [WIP] Sincei core library

Backend library for [sincei](https://github.com/bhardwaj-lab/sincei).

Temporarily includes a [typer](https://github.com/fastapi/typer) CLI to run commands.

To install in current environment use [maturin](https://www.maturin.rs/installation.html).
First, create a new conda environment, then clone the git repository, and lastly compile the
program using maturin.

```bash
conda create -n sincei_rust -c conda-forge python=3.12 maturin

git clone https://github.com/bhardwaj-lab/sincei_rust.git

maturin build --release
```
