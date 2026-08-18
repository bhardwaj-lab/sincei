#!/usr/bin/env python
"""Snapshot tests for ``scFilterStats`` (rust sincei).

Each flag/option is exercised on the data in ``tests/testdata/`` and the tool's
per-cell TSV output is compared against a committed snapshot under
``tests/testdata/scFilterStats/``.

Most scenarios use dense sampling (``-bs 2000 --distance-between-bins 0``) so the
tiny test BAM's reads are actually captured; the default 100kb-every-1.1Mb grid
steps over them and yields all-zero stats.  The ``bin_size`` and
``distance_between_bins`` scenarios set their own sampling (they *are* the
sampling flags).

Run the tests::

    pytest tests/test_scFilterStats.py

Regenerate the snapshots after an intentional output change::

    python tests/test_scFilterStats.py --update
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "testdata"
SNAP_DIR = DATA / "scFilterStats"
BAM = str(DATA / "test_i1.bam")
BARCODES = str(DATA / "test_barcodes.txt")
BED = str(DATA / "Chrna9_regions.bed")

TOOL = "scFilterStats"
BASE = ["-b", BAM, "-bc", BARCODES, "-ct", "BC", "-p", "1"]

# Dense sampling so the reads are captured; prepended to every scenario that
# does not set its own sampling.
SAMPLING = ["-bs", "2000", "--distanceBetweenBins", "0"]

# name -> (extra args, own_sampling).  When own_sampling is False, SAMPLING is
# prepended.  Every comparable flag is covered; --motifFilter/--genome2bit are
# omitted (no 2bit genome in testdata).
SCENARIOS: dict[str, tuple[list[str], bool]] = {
    "baseline": ([], False),
    "bin_size": (["-bs", "5000", "--distanceBetweenBins", "0"], True),
    "distance_between_bins": (["-bs", "2000", "--distanceBetweenBins", "1000"], True),
    "dup_start_bc": (["--duplicateFilter", "start_bc"], False),
    "dup_start_end_bc": (["--duplicateFilter", "start_end_bc"], False),
    "dup_umi": (["--duplicateFilter", "start_bc_umi"], False),
    "gc_content": (["-gc", "0.2,0.8"], False),
    "min_aligned_fraction": (["--minAlignedFraction", "0.5"], False),
    "min_mapping_quality": (["-mq", "20"], False),
    "sam_flag_include": (["--samFlagInclude", "64"], False),
    "sam_flag_exclude": (["--samFlagExclude", "16"], False),
    "filter_rna_strand": (["--filterRNAstrand", "forward"], False),
    "labels": (["-l", "mysample"], False),
    "smart_labels": (["--smartLabels"], False),
    "blacklist": (["-bl", BED], False),
}


def _args(name: str) -> list[str]:
    extra, own = SCENARIOS[name]
    return extra if own else SAMPLING + extra


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


# The exact output contract. Spelled out here rather than derived from the code,
# so that changing the tool's columns fails this test instead of silently
# rewriting every snapshot.
EXPECTED_HEADER = (
    "Cell_ID\tTotal_sampled\tFiltered\tBlacklisted\tLow_MAPQ\tMissing_Flags\tExcluded_Flags\t"
    "Internal_Duplicates\tMarked_Duplicates\tSingletons\tWrong_strand\tWrong_motif\t"
    "Unwanted_GC_content\tLow_aligned_fraction"
)


def test_header_is_exact(tmp_path: Path) -> None:
    text = run_tool(_args("baseline"), str(tmp_path / "header.tsv"))
    assert text.splitlines()[0] == EXPECTED_HEADER


def test_cell_id_joins_sample_and_barcode(tmp_path: Path) -> None:
    """Cell_ID is `{sample}::{barcode}`, one column, not a sample/barcode pair."""
    text = run_tool(["-l", "mysample", *SAMPLING], str(tmp_path / "cell_id.tsv"))
    header, *rows = text.splitlines()

    assert header.split("\t")[0] == "Cell_ID"
    barcodes = [
        line.strip() for line in Path(BARCODES).read_text().splitlines() if line.strip()
    ]
    cell_ids = [row.split("\t")[0] for row in rows]
    assert cell_ids == [f"mysample::{bc}" for bc in barcodes]

    # One Cell_ID column + 13 stat columns, on the header and every row alike.
    assert len(header.split("\t")) == 14
    assert all(len(row.split("\t")) == 14 for row in rows)


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_snapshot(name: str, tmp_path: Path) -> None:
    snap = SNAP_DIR / f"{name}.tsv"
    if not snap.exists():
        pytest.skip(
            f"missing snapshot {snap} (regenerate: python {Path(__file__).name} --update)"
        )
    got = run_tool(_args(name), str(tmp_path / f"{name}.tsv"))
    assert _norm(got) == _norm(snap.read_text()), (
        f"output for '{name}' differs from snapshot"
    )


def _generate() -> None:
    SNAP_DIR.mkdir(parents=True, exist_ok=True)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        for name in SCENARIOS:
            text = run_tool(_args(name), f"{tmp}/{name}.tsv")
            (SNAP_DIR / f"{name}.tsv").write_text(text)
            print(f"wrote {SNAP_DIR / f'{name}.tsv'}")


if __name__ == "__main__":
    if "--update" in sys.argv:
        _generate()
    else:
        raise SystemExit("use --update to regenerate snapshots, or run via pytest")
