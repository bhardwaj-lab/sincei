#!/usr/bin/env python
"""Snapshot tests for ``scCountReads`` (rust sincei).

Each scenario runs the tool on ``tests/testdata/`` and compares a rendering of
the resulting AnnData against a committed snapshot under
``tests/testdata/scCountReads/``.  The h5ad itself is not compared: HDF5 is a
binary container whose bytes vary between writes, so :func:`render` projects it
to the part that carries meaning -- the shape, and every non-zero count keyed by
cell and region.

``scCountReads`` has two subcommands and they share almost every flag, so the
scenarios are declared once and the ones that apply to both are run twice.

Not covered: ``--motifFilter`` / ``--genome2bit`` need a 2bit genome that is not
in ``testdata/``.

Run the tests::

    pytest tests/test_scCountReads.py
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import anndata as ad
import pytest
from _cli_testing import (
    BAM1,
    BAM2,
    BAM_MERGED,
    BARCODES,
    BED,
    DATA,
    Scenario,
    norm,
    read_snapshot,
    run,
    tool_path,
)

if TYPE_CHECKING:
    from pathlib import Path

TOOL = "scCountReads"
GTF = str(DATA / "Chrna9.gtf")

# One BAM keeps the snapshots small; the multi-BAM and merged cases have their
# own tests below.
BASE = ["-b", BAM1, "-bc", BARCODES, "-ct", "BC", "-p", "1"]
BINS = ["-bs", "100000"]
FEATURES = ["--bed", BED]


def render(path: str) -> str:
    """Project an AnnData to the text a snapshot can meaningfully diff.

    Shape first, then one line per non-zero count as ``cell region value``,
    sorted so the rendering does not depend on storage order.
    """
    adata = ad.read_h5ad(path)
    rows, cols = adata.shape
    lines = [f"shape\t{rows}\t{cols}"]

    entries = [
        (str(adata.obs_names[i]), str(adata.var_names[j]), int(v))
        for i, j, v in zip(*adata.X.nonzero(), adata.X.data, strict=False)
    ]
    lines += [f"{cell}\t{region}\t{value}" for cell, region, value in sorted(entries)]
    return "\n".join(lines)


# name -> scenario.  `own_sampling` marks the ones that set their own --binSize,
# so the default is not prepended.
SCENARIOS: dict[str, Scenario] = {
    "baseline": Scenario(
        [], differs=False, reason="it is the reference the others are compared against"
    ),
    "bin_size": Scenario(["-bs", "20000"], own_sampling=True),
    "region": Scenario(["-r", "5:65900000-66000000"]),
    "chr_to_skip": Scenario(["--chrToSkip", "5"]),
    "min_mapping_quality": Scenario(["-mq", "40"]),
    "blacklist": Scenario(["-bl", BED]),
    "sam_flag_exclude": Scenario(["--samFlagExclude", "16"]),
    "gc_content": Scenario(["-gc", "0.30,0.50"]),
    "min_aligned_fraction": Scenario(["--minAlignedFraction", "0.8"]),
    "filter_rna_strand": Scenario(["--filterRNAstrand", "forward"]),
    "duplicate_filter": Scenario(["--duplicateFilter", "start_bc"]),
    # Read adjustments only change which bin a read lands in when the bins
    # are comparable to the adjustment, so these set their own fine binning.
    "extend_reads": Scenario(["-bs", "100", "-e", "300"], own_sampling=True),
    # Centring a read on its own span is a no-op, so it is only observable
    # once the span has been extended.
    "center_reads": Scenario(
        ["-bs", "100", "-e", "300", "--centerReads"], own_sampling=True
    ),
    "labels": Scenario(["-l", "MYSAMPLE"]),
}


def _args(name: str, mode: str) -> list[str]:
    scenario = SCENARIOS[name]
    if mode == "features":
        mode_args = FEATURES
    else:
        mode_args = [] if scenario.own_sampling else BINS
    return [*mode_args, *scenario.args]


def _count(mode: str, base: list[str], extra: list[str], out: Path) -> str:
    """Run a subcommand, assert it succeeded, return the rendered AnnData.

    Not :func:`run_ok`: that reads the output file as text to return it, which
    an h5ad is not.
    """
    proc = run(tool_path(TOOL), [mode, *base, "-o", str(out), *extra])
    assert proc.returncode == 0, (
        f"{TOOL} {mode} {' '.join(extra)} failed "
        f"(exit {proc.returncode}):\n{proc.stderr}"
    )
    return render(str(out))


def _run(mode: str, extra: list[str], out: Path) -> str:
    return _count(mode, BASE, extra, out)


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_bins_snapshot(name: str, tmp_path: Path) -> None:
    got = _run("bins", _args(name, "bins"), tmp_path / "out.h5ad")
    assert norm(got) == norm(read_snapshot(TOOL, f"bins_{name}")), (
        f"bins output for '{name}' differs from its snapshot"
    )


@pytest.mark.parametrize("name", ["baseline", "region", "min_mapping_quality"])
def test_features_snapshot(name: str, tmp_path: Path) -> None:
    got = _run("features", _args(name, "features"), tmp_path / "out.h5ad")
    assert norm(got) == norm(read_snapshot(TOOL, f"features_{name}")), (
        f"features output for '{name}' differs from its snapshot"
    )


def test_every_scenario_that_should_differ_does() -> None:
    """A scenario whose snapshot equals baseline asserts nothing about its flag."""
    baseline = norm(read_snapshot(TOOL, "bins_baseline"))
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            continue
        assert norm(read_snapshot(TOOL, f"bins_{name}")) != baseline, (
            f"scenario '{name}' ({' '.join(scenario.args)}) produces the baseline "
            f"output, so it cannot fail if the flag breaks. Pick a discriminating "
            f"value, or mark it differs=False with a reason."
        )


def test_scenarios_that_cannot_differ_say_why() -> None:
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            assert scenario.reason, f"scenario '{name}' needs a recorded reason"


# The two subcommands


def test_features_counts_only_the_annotated_regions(tmp_path: Path) -> None:
    """``bins`` tiles the whole chromosome; ``features`` uses the BED alone."""
    bins = _run("bins", [*BINS], tmp_path / "bins.h5ad")
    features = _run("features", [*FEATURES], tmp_path / "features.h5ad")

    bins_cols = int(bins.splitlines()[0].split("\t")[2])
    feature_cols = int(features.splitlines()[0].split("\t")[2])
    assert feature_cols < bins_cols


def test_metagene_groups_exons_onto_their_gene(tmp_path: Path) -> None:
    per_transcript = _run("features", ["--bed", GTF], tmp_path / "tx.h5ad")
    per_gene = _run("features", ["--bed", GTF, "--metagene"], tmp_path / "gene.h5ad")

    cols = lambda text: int(text.splitlines()[0].split("\t")[2])  # noqa: E731
    assert cols(per_gene) <= cols(per_transcript)


# Sample labelling


def test_labels_name_the_cells_and_default_to_the_file_stem(tmp_path: Path) -> None:
    default = _run("bins", [*BINS], tmp_path / "default.h5ad")
    labelled = _run("bins", [*BINS, "-l", "MYSAMPLE"], tmp_path / "labelled.h5ad")

    assert "test_i1::" in default
    assert "MYSAMPLE::" in labelled
    assert "test_i1::" not in labelled


def test_a_label_per_bam_is_required(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        [
            "bins",
            "-b",
            BAM1,
            BAM2,
            "-bc",
            BARCODES,
            "-ct",
            "BC",
            "-o",
            str(tmp_path / "out.h5ad"),
            *BINS,
            "-l",
            "onlyone",
        ],
    )
    assert proc.returncode != 0
    assert "--labels count" in proc.stdout + proc.stderr


# Merged BAMs: --groupTag


def test_a_merged_bam_grouped_by_tag_matches_its_separate_sources(
    tmp_path: Path,
) -> None:
    """The invariant --groupTag exists to satisfy: one merged file split back
    apart by read group must equal counting its two sources separately."""
    separate = _count(
        "bins",
        ["-b", BAM1, BAM2, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*BINS],
        tmp_path / "separate.h5ad",
    )
    grouped = _count(
        "bins",
        ["-b", BAM_MERGED, "-bc", BARCODES, "-ct", "BC", "-p", "1"],
        [*BINS, "-gt", "RG"],
        tmp_path / "grouped.h5ad",
    )

    assert grouped == separate


def test_a_barcode_shared_across_samples_stays_two_cells(tmp_path: Path) -> None:
    """``GCGAGCAT`` is in both sources, so per-barcode counting merges two cells."""
    merged = ["-b", BAM_MERGED, "-bc", BARCODES, "-ct", "BC", "-p", "1"]

    flat = _count("bins", merged, [*BINS], tmp_path / "flat.h5ad")
    grouped = _count("bins", merged, [*BINS, "-gt", "RG"], tmp_path / "grouped.h5ad")

    def shared(text: str) -> set[str]:
        return {
            line.split("\t")[0]
            for line in text.splitlines()[1:]
            if line.split("\t")[0].endswith("GCGAGCAT")
        }

    assert shared(flat) == {"test_i1_i2::GCGAGCAT"}
    assert shared(grouped) == {
        "test_i1::GCGAGCAT",
        "test_i2::GCGAGCAT",
    }


def test_group_tag_with_several_bams_is_rejected(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        [
            "bins",
            "-b",
            BAM1,
            BAM2,
            "-bc",
            BARCODES,
            "-ct",
            "BC",
            "-o",
            str(tmp_path / "out.h5ad"),
            *BINS,
            "-gt",
            "RG",
        ],
    )
    assert proc.returncode != 0
    assert "single merged BAM" in proc.stdout + proc.stderr


def test_group_tag_on_a_bam_without_read_groups_is_rejected(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        ["bins", *BASE, "-o", str(tmp_path / "out.h5ad"), *BINS, "-gt", "RG"],
    )
    assert proc.returncode != 0
    assert "@RG" in proc.stdout + proc.stderr


# Error paths


def test_a_missing_bam_fails_with_a_message_naming_it(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        [
            "bins",
            "-b",
            "/nonexistent/reads.bam",
            "-bc",
            BARCODES,
            "-ct",
            "BC",
            "-o",
            str(tmp_path / "out.h5ad"),
            *BINS,
        ],
    )
    assert proc.returncode != 0
    assert "/nonexistent/reads.bam" in proc.stdout + proc.stderr


def test_a_barcode_tag_the_bam_lacks_fails_with_advice(tmp_path: Path) -> None:
    proc = run(
        tool_path(TOOL),
        [
            "bins",
            "-b",
            BAM1,
            "-bc",
            BARCODES,
            "-ct",
            "ZZ",
            "-o",
            str(tmp_path / "out.h5ad"),
            *BINS,
        ],
    )
    assert proc.returncode != 0
    assert "ZZ" in proc.stdout + proc.stderr
