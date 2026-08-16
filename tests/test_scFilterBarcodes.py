#!/usr/bin/env python
"""Snapshot tests for ``scFilterBarcodes`` (rust sincei).

Each flag/option is exercised on the data in ``tests/testdata/`` and the tool's
TSV output is compared against a committed snapshot under
``tests/testdata/scFilterBarcodes/``.  The snapshots were validated to be
equivalent to the original Python ``sincei`` implementation.

The ``-rp`` (rank plot) scenario is special: it only checks that the tool runs
and writes a PNG (the image itself is not snapshotted).

Run the tests::

    pytest tests/test_scFilterBarcodes.py

Regenerate the snapshots after an intentional output change::

    python tests/test_scFilterBarcodes.py --update
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "testdata"
SNAP_DIR = DATA / "scFilterBarcodes"
BAM = str(DATA / "test_i1.bam")
WHITELIST = str(DATA / "test_barcodes.txt")
BED = str(DATA / "Chrna9_regions.bed")

TOOL = "scFilterBarcodes"
BASE = ["-b", BAM, "-w", WHITELIST, "-ct", "BC", "-p", "1"]

# name -> extra args.  Short flags are the tool's own; every flag is covered:
# -b/-w/-o/-ct/-p via BASE, plus these.  (-rp is tested separately.)
SCENARIOS: dict[str, list[str]] = {
    "baseline": [],
    "min_count": ["-mc", "2"],
    "min_hamming_dist": ["-d", "1"],
    "min_mapping_quality": ["-mq", "20"],
    "blacklist": ["-bl", BED],
    "bin_size": ["-bs", "50000"],
    "group_tag": ["-gt", "RG"],
    "verbose": ["-v"],
}


def _tool_path() -> str:
    path = shutil.which(TOOL)
    if path is None:
        pytest.skip(f"{TOOL} not on PATH (run pytest inside the sincei_dev env)")
    return path


def run_tool(extra: list[str], out: str) -> str:
    """Run the tool, assert it succeeds, return the output TSV text."""
    cmd = [_tool_path(), *BASE, "-o", out, *extra]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"{TOOL} {extra} failed:\n{proc.stderr}"
    return Path(out).read_text()


def _norm(text: str) -> str:
    """Normalize for comparison: strip trailing whitespace, drop blank tail."""
    return "\n".join(line.rstrip() for line in text.splitlines()).rstrip("\n")


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_snapshot(name: str, tmp_path: Path) -> None:
    snap = SNAP_DIR / f"{name}.tsv"
    if not snap.exists():
        pytest.skip(
            f"missing snapshot {snap} (regenerate: python {Path(__file__).name} --update)"
        )
    got = run_tool(SCENARIOS[name], str(tmp_path / f"{name}.tsv"))
    assert _norm(got) == _norm(snap.read_text()), (
        f"output for '{name}' differs from snapshot"
    )


def test_rank_plot_writes_png(tmp_path: Path) -> None:
    """`-rp` should run and write a PNG; the image is not snapshotted."""
    png = tmp_path / "rank.png"
    cmd = [_tool_path(), *BASE, "-o", str(tmp_path / "out.tsv"), "-rp", str(png)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"{TOOL} -rp failed:\n{proc.stderr}"
    assert png.exists() and png.stat().st_size > 0, "rank plot PNG was not written"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "output is not a valid PNG"


def _generate() -> None:
    SNAP_DIR.mkdir(parents=True, exist_ok=True)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        for name, extra in SCENARIOS.items():
            text = run_tool(extra, f"{tmp}/{name}.tsv")
            (SNAP_DIR / f"{name}.tsv").write_text(text)
            print(f"wrote {SNAP_DIR / f'{name}.tsv'}")


if __name__ == "__main__":
    if "--update" in sys.argv:
        _generate()
    else:
        raise SystemExit("use --update to regenerate snapshots, or run via pytest")
