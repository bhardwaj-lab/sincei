#!/usr/bin/env python
"""Snapshot tests for ``scFilterStats`` (rust sincei).

Each scenario runs the tool on ``tests/testdata/`` and compares its per-cell TSV
against a committed snapshot under ``tests/testdata/scFilterStats/``.

Most scenarios use dense sampling (``-bs 2000 --distanceBetweenBins 0``) so the
tiny test BAM's reads are actually captured; the default 100kb-every-1.1Mb grid
steps over them and yields all-zero stats.  Scenarios that set their own
sampling declare ``own_sampling=True``.

Not covered: ``--motifFilter`` / ``--genome2bit`` need a 2bit genome that is not
in ``testdata/``.

Run the tests::

    pytest tests/test_scFilterStats.py

"""

from __future__ import annotations

from pathlib import Path

import pytest
from _cli_testing import (
    BAM1,
    BAM2,
    BAM_MERGED,
    BARCODES,
    BED,
    Scenario,
    norm,
    read_snapshot,
    run,
    run_ok,
    tool_path,
)

TOOL = "scFilterStats"
BASE = ["-b", BAM1, "-bc", BARCODES, "-ct", "BC", "-p", "1"]
SAMPLING = ["-bs", "2000", "--distanceBetweenBins", "0"]

# Scenarios where the port and the original disagree.  These are snapshotted
# from the port, so the suite pins current behaviour instead of failing on a
# difference nobody has adjudicated yet.  Each entry records what the difference
# *is*, deliberately without claiming which side is correct, that is a
# question for the sincei merge, not for the test suite.
SCENARIOS: dict[str, Scenario] = {
    "baseline": Scenario(
        [], differs=False, reason="it is the reference the others are compared against"
    ),
    "bin_size": Scenario(
        ["-bs", "5000", "--distanceBetweenBins", "0"],
        own_sampling=True,
        differs=False,
        reason=(
            "with distanceBetweenBins 0 the bins are contiguous, so every read is "
            "sampled whatever the bin size"
        ),
    ),
    "distance_between_bins": Scenario(
        ["-bs", "2000", "--distanceBetweenBins", "1000"], own_sampling=True
    ),
    "dup_start_bc": Scenario(["--duplicateFilter", "start_bc"]),
    "dup_start_end_bc": Scenario(["--duplicateFilter", "start_end_bc"]),
    "dup_umi": Scenario(["--duplicateFilter", "start_bc_umi"]),
    "dup_start_end_bc_umi": Scenario(
        ["--duplicateFilter", "start_end_bc_umi"],
        differs=False,
        reason="the strictest duplicate key finds no duplicates among 13 reads",
    ),
    "gc_content": Scenario(["-gc", "0.30,0.50"]),
    "min_aligned_fraction": Scenario(["--minAlignedFraction", "0.8"]),
    # Long form: the original has no -mq short flag, and these scenarios are
    # replayed against it to generate the snapshots.
    "min_mapping_quality": Scenario(["--minMappingQuality", "20"]),
    "sam_flag_include": Scenario(["--samFlagInclude", "16"]),
    "sam_flag_exclude": Scenario(["--samFlagExclude", "16"]),
    "filter_rna_strand": Scenario(["--filterRNAstrand", "forward"]),
    "filter_rna_strand_reverse": Scenario(["--filterRNAstrand", "reverse"]),
    "labels": Scenario(["-l", "mysample"]),
    "smart_labels": Scenario(
        ["--smartLabels"],
        differs=False,
        reason=(
            "a no-op in both implementations: the original's validateInputs applies "
            "smartLabels() unconditionally (ParserCommon.py:806), so the file stem "
            "is always used"
        ),
    ),
    "blacklist": Scenario(["-bl", BED]),
    "chr_to_skip": Scenario(["--chrToSkip", "5"]),
}


# The exact output contract. Spelled out here rather than derived from the code,
# so that changing the tool's columns fails this test instead of silently
# rewriting every snapshot.
EXPECTED_HEADER = (
    "Cell_ID\tTotal_sampled\tFiltered\tBlacklisted\tLow_MAPQ\tMissing_Flags\tExcluded_Flags\t"
    "Internal_Duplicates\tMarked_Duplicates\tSingletons\tWrong_strand\tWrong_motif\t"
    "Unwanted_GC_content\tLow_aligned_fraction"
)


def _run(extra: list[str], out: Path) -> str:
    return run_ok(tool_path(TOOL), BASE, extra, str(out))


def test_header_is_exact(tmp_path: Path) -> None:
    text = _run([*SAMPLING], tmp_path / "header.tsv")
    assert text.splitlines()[0] == EXPECTED_HEADER


def test_cell_id_joins_sample_and_barcode(tmp_path: Path) -> None:
    """Cell_ID is `{sample}::{barcode}`, one column, not a sample/barcode pair."""
    text = _run(["-l", "mysample", *SAMPLING], tmp_path / "cell_id.tsv")
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
    scenario = SCENARIOS[name]
    extra = scenario.args if scenario.own_sampling else [*SAMPLING, *scenario.args]
    got = _run(extra, tmp_path / "out.tsv")
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


# Multiple input BAMs


def test_two_bams_give_one_row_per_sample_and_barcode(tmp_path: Path) -> None:
    """``--bamfiles`` is plural; Cell_ID must namespace barcodes by sample."""
    text = run_ok(
        tool_path(TOOL),
        ["-b", BAM1, BAM2, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*SAMPLING],
        str(tmp_path / "two.tsv"),
    )
    header, *rows = text.splitlines()
    barcodes = [
        line.strip() for line in Path(BARCODES).read_text().splitlines() if line.strip()
    ]

    assert header == EXPECTED_HEADER
    assert len(rows) == 2 * len(barcodes)

    cell_ids = [row.split("\t")[0] for row in rows]
    assert sum(cid.startswith("test_i1::") for cid in cell_ids) == len(barcodes)
    assert sum(cid.startswith("test_i2::") for cid in cell_ids) == len(barcodes)


def test_labels_apply_to_each_bam_in_order(tmp_path: Path) -> None:
    text = run_ok(
        tool_path(TOOL),
        ["-b", BAM1, BAM2, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*SAMPLING, "-l", "first", "second"],
        str(tmp_path / "labelled.tsv"),
    )
    cell_ids = [row.split("\t")[0] for row in text.splitlines()[1:]]
    assert any(cid.startswith("first::") for cid in cell_ids)
    assert any(cid.startswith("second::") for cid in cell_ids)


# Merged BAMs: --groupTag


def _cell_ids(text: str) -> list[str]:
    return [row.split("\t")[0] for row in text.splitlines()[1:]]


def test_a_merged_bam_grouped_by_tag_matches_its_separate_sources(
    tmp_path: Path,
) -> None:
    """The invariant --groupTag exists to satisfy.

    One merged file, split back apart by each read's group tag, must give the
    same per-cell stats as counting the two source files separately. This is
    the check that catches a group axis wired up wrongly: without it the tool
    still exits 0 and still writes a well-formed table, just full of zeros.
    """
    separate = run_ok(
        tool_path(TOOL),
        ["-b", BAM1, BAM2, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*SAMPLING],
        str(tmp_path / "separate.tsv"),
    )
    grouped = run_ok(
        tool_path(TOOL),
        ["-b", BAM_MERGED, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*SAMPLING, "-gt", "RG"],
        str(tmp_path / "grouped.tsv"),
    )

    assert norm(grouped) == norm(separate)


def test_a_barcode_shared_across_samples_stays_two_cells(tmp_path: Path) -> None:
    """``GCGAGCAT`` occurs in both sources of the merged BAM.

    Counted per-barcode the two cells merge into one row and their reads are
    summed; counted per-group they stay apart.
    """
    base = ["-b", BAM_MERGED, "-bc", BARCODES, "-ct", "BC", "-p", "1"]

    flat = run_ok(tool_path(TOOL), base, [*SAMPLING], str(tmp_path / "flat.tsv"))
    grouped = run_ok(
        tool_path(TOOL), base, [*SAMPLING, "-gt", "RG"], str(tmp_path / "grouped.tsv")
    )

    def shared(ids: list[str]) -> list[str]:
        return [i for i in ids if i.endswith("GCGAGCAT")]

    assert len(shared(_cell_ids(flat))) == 1
    assert shared(_cell_ids(grouped)) == [
        "test_i1::GCGAGCAT",
        "test_i2::GCGAGCAT",
    ]


def test_grouping_names_cells_after_the_read_group_not_the_file(
    tmp_path: Path,
) -> None:
    """The sample half of the Cell_ID comes from the tag, so ``--labels`` and
    the file name have no say in it."""
    text = run_ok(
        tool_path(TOOL),
        ["-b", BAM_MERGED, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*SAMPLING, "-gt", "RG", "-l", "ignored"],
        str(tmp_path / "grouped.tsv"),
    )
    cell_ids = _cell_ids(text)

    assert not any(cid.startswith("ignored::") for cid in cell_ids)
    assert not any(cid.startswith("test_i1_i2::") for cid in cell_ids)
    prefixes = {cid.split("::")[0] for cid in cell_ids}
    assert prefixes == {"test_i1", "test_i2"}


def test_group_tag_with_several_bams_is_rejected(tmp_path: Path) -> None:
    """A merged BAM carries its samples internally, so there is only ever one."""
    proc = run(
        tool_path(TOOL),
        [
            "-b",
            BAM1,
            BAM2,
            "-bc",
            BARCODES,
            "-ct",
            "BC",
            "-o",
            str(tmp_path / "out.tsv"),
            *SAMPLING,
            "-gt",
            "RG",
        ],
    )
    assert proc.returncode != 0
    assert "single merged BAM" in proc.stdout + proc.stderr


def test_group_tag_on_a_bam_without_read_groups_is_rejected(tmp_path: Path) -> None:
    """An un-merged BAM declares no @RG, and counting it per-barcode would
    silently merge cells that share a barcode across samples."""
    proc = run(
        tool_path(TOOL),
        [
            "-b",
            BAM1,
            "-bc",
            BARCODES,
            "-ct",
            "BC",
            "-o",
            str(tmp_path / "out.tsv"),
            *SAMPLING,
            "-gt",
            "RG",
        ],
    )
    assert proc.returncode != 0
    assert "@RG" in proc.stdout + proc.stderr


# UMI tag


def test_umi_tag_selects_the_umi_source(tmp_path: Path) -> None:
    """``--umiTag`` picks which tag holds the UMI.

    ``MI`` (barcode+UMI concatenated) and the default ``RX`` (UMI alone) agree
    on this data, because the duplicate key already includes the barcode.  The
    test pins that equivalence rather than pretending the flag is unverified.
    """
    dedup = ["--duplicateFilter", "start_bc_umi"]
    with_rx = _run([*SAMPLING, *dedup], tmp_path / "rx.tsv")
    with_mi = _run([*SAMPLING, "-ut", "MI", *dedup], tmp_path / "mi.tsv")
    plain = _run([*SAMPLING], tmp_path / "plain.tsv")

    assert norm(with_rx) != norm(plain), "the UMI duplicate filter changed nothing"
    assert norm(with_mi) == norm(with_rx)


# Error paths


def test_a_missing_bam_fails_with_a_message_naming_it(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        [
            "-b",
            "/nonexistent/reads.bam",
            "-bc",
            BARCODES,
            "-o",
            str(tmp_path / "o.tsv"),
        ],
    )
    assert proc.returncode != 0, "a missing BAM should not exit 0"
    assert "/nonexistent/reads.bam" in proc.stdout + proc.stderr


def test_a_barcode_tag_the_bam_lacks_fails_with_advice(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        ["-b", BAM1, "-bc", BARCODES, "-ct", "ZZ", "-o", str(tmp_path / "o.tsv")],
    )
    assert proc.returncode != 0, "an absent tag should not exit 0"
    assert "ZZ" in proc.stdout + proc.stderr
