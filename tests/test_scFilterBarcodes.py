#!/usr/bin/env python
"""Snapshot tests for ``scFilterBarcodes`` (rust sincei).

Each scenario runs the tool on ``tests/testdata/`` and compares its TSV output
against a committed snapshot under ``tests/testdata/scFilterBarcodes/``.

Known limits of this fixture, so the gaps are not mistaken for coverage:

* The tool counts *bins occupied per barcode*, and every barcode's reads in
  ``test_i1.bam`` are co-located (mate pairs share a start).  So every count is
  1, which pins ``count_log10`` at 0.0 and ``count_rank`` at 1 in every
  scenario.  Neither column, nor ``--minCount`` as anything but an
  all-or-nothing threshold, is verified.  ``--binSize`` cannot change the
  output for the same reason.
* ``--motifFilter`` / ``--genome2bit`` need a 2bit genome that is not in
  ``testdata/``.

Run the tests::

    pytest tests/test_scFilterBarcodes.py

"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from _cli_testing import (
    BAM1,
    BAM_MERGED,
    BARCODES,
    BARCODES_1MIS,
    BED,
    Scenario,
    norm,
    read_snapshot,
    run,
    run_ok,
    tool_path,
)

if TYPE_CHECKING:
    from pathlib import Path

TOOL = "scFilterBarcodes"
BASE = ["-b", BAM1, "-w", BARCODES, "-ct", "BC", "-p", "1"]

SCENARIOS: dict[str, Scenario] = {
    "baseline": Scenario(
        [], differs=False, reason="it is the reference the others are compared against"
    ),
    "min_count": Scenario(["-mc", "2"]),
    "min_mapping_quality": Scenario(["-mq", "40"]),
    "blacklist": Scenario(["-bl", BED]),
    "chr_to_skip": Scenario(["--chrToSkip", "5"]),
    "bin_size": Scenario(
        ["-bs", "50000"],
        differs=False,
        reason=(
            "every barcode's reads are co-located, so each occupies exactly one "
            "bin at any bin size (see the module docstring)"
        ),
    ),
    "verbose": Scenario(
        ["-v"], differs=False, reason="affects stderr only, never the TSV"
    ),
}


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_snapshot(name: str, tmp_path: Path) -> None:
    got = run_ok(tool_path(TOOL), BASE, SCENARIOS[name].args, str(tmp_path / "out.tsv"))
    assert norm(got) == norm(read_snapshot(TOOL, name)), (
        f"output for '{name}' differs from its snapshot"
    )


def test_every_scenario_that_should_differ_does() -> None:
    """A scenario whose snapshot equals baseline asserts nothing about its flag.

    This guards the suite against the failure it was written to fix: a value
    that looks like it exercises a flag but does not discriminate on this data.
    """
    baseline = norm(read_snapshot(TOOL, "baseline"))
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            continue
        assert norm(read_snapshot(TOOL, name)) != baseline, (
            f"scenario '{name}' ({' '.join(scenario.args)}) produces the baseline "
            f"output, so it cannot fail if the flag breaks. Pick a discriminating "
            f"value, or mark it differs=False with a reason."
        )


def test_scenarios_that_cannot_differ_say_why() -> None:
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            assert scenario.reason, f"scenario '{name}' needs a recorded reason"


# Fuzzy whitelist matching


def test_exact_matching_rejects_a_one_mismatch_whitelist(tmp_path: Path) -> None:
    """``test_barcodes_1mis.txt`` is every barcode with one base substituted."""
    text = run_ok(
        tool_path(TOOL),
        ["-b", BAM1, "-w", BARCODES_1MIS, "-ct", "BC", "-p", "1"],
        ["-d", "0"],
        str(tmp_path / "exact.tsv"),
    )
    assert len(text.splitlines()) == 1, "only the header should remain"


def test_fuzzy_matching_recovers_a_one_mismatch_whitelist(tmp_path: Path) -> None:
    text = run_ok(
        tool_path(TOOL),
        ["-b", BAM1, "-w", BARCODES_1MIS, "-ct", "BC", "-p", "1"],
        ["-d", "1"],
        str(tmp_path / "fuzzy.tsv"),
    )
    rows = text.splitlines()[1:]
    assert len(rows) == 8, f"expected the 8 barcodes of {BAM1}, got {len(rows)}"


# Rank plot


def test_rank_plot_writes_png_even_when_every_count_is_equal(tmp_path: Path) -> None:
    """``-rp`` must survive uniform counts.

    This is exactly the input that crashes the original sincei: all counts are
    1, so ``count_log10`` is uniformly 0 and its axis construction raises
    ``ValueError: arange: cannot compute length``.
    """
    png = tmp_path / "rank.png"
    proc = run(
        tool_path(TOOL),
        [*BASE, "-o", str(tmp_path / "out.tsv"), "-rp", str(png)],
    )
    assert proc.returncode == 0, f"{TOOL} -rp failed:\n{proc.stderr}"
    assert png.exists(), "rank plot PNG was not written"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "output is not a valid PNG"


# Merged BAMs: --groupTag


def test_a_shared_barcode_becomes_two_entries_under_group_tag(tmp_path: Path) -> None:
    """``GCGAGCAT`` is in both source samples of the merged BAM.

    Without the flag the two cells collapse into one entry; with it they stay
    apart, which is the whole reason the option exists.
    """
    base = ["-b", BAM_MERGED, "-w", BARCODES, "-ct", "BC", "-p", "1"]

    flat = run_ok(tool_path(TOOL), base, [], str(tmp_path / "flat.tsv"))
    grouped = run_ok(
        tool_path(TOOL), base, ["-gt", "RG"], str(tmp_path / "grouped.tsv")
    )

    def names(text: str) -> list[str]:
        return [line.split("\t")[0] for line in text.splitlines()[1:]]

    assert [n for n in names(flat) if n.endswith("GCGAGCAT")] == ["GCGAGCAT"]
    assert [n for n in names(grouped) if n.endswith("GCGAGCAT")] == [
        "test_i1::GCGAGCAT",
        "test_i2::GCGAGCAT",
    ]
    assert len(names(grouped)) == len(names(flat)) + 1


def test_group_tag_on_a_bam_without_read_groups_is_rejected(tmp_path: Path) -> None:
    """An un-merged BAM has no @RG, and counting it per-barcode would be wrong."""
    proc = run(
        tool_path(TOOL),
        [*BASE, "-o", str(tmp_path / "out.tsv"), "-gt", "RG"],
    )
    assert proc.returncode != 0
    assert "@RG" in proc.stdout + proc.stderr


# Error paths


def test_a_missing_bam_fails_with_a_message_naming_it(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        ["-b", "/nonexistent/reads.bam", "-ct", "BC", "-o", str(tmp_path / "o.tsv")],
    )
    assert proc.returncode != 0, "a missing BAM should not exit 0"
    assert "/nonexistent/reads.bam" in proc.stdout + proc.stderr


def test_a_barcode_tag_the_bam_lacks_fails_with_advice(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        ["-b", BAM1, "-ct", "ZZ", "-o", str(tmp_path / "o.tsv")],
    )
    assert proc.returncode != 0, "an absent tag should not exit 0"
    assert "ZZ" in proc.stdout + proc.stderr
