# [WIP] Sincei core library

Backend library for [sincei](https://github.com/bhardwaj-lab/sincei).

Temporarily includes a [typer](https://github.com/fastapi/typer) CLI to run commands.

To install in current environment we need _rust_, _python_ and _maturin_.

You can install rust via [rustup](http://rust-lang.org/tools/install/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Download _python_ and _maturin_ from conda. First, create a new conda environment, then
clone this git repository, and lastly compile the program using maturin.

```bash
conda create -n sincei_rust -c conda-forge python=3.12 maturin cmake
git clone https://github.com/bhardwaj-lab/sincei_rust.git
maturin build --release
pip install -e .
```

To update the extension, pull the changes from the repo and recompile it:
```bash
giut pull
maturin build --release
```
