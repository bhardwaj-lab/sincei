"""Snapshot tests for ``scBulkCoverage`` (rust sincei).

Each scenario runs the tool on the SL2 fixture (two paired-end BAMs, ten cells
in two clusters, the Ogfrl1 locus on chr1) and compares every output bigwig,
rendered to text by :func:`render`, with a committed snapshot under
``tests/testdata/scBulkCoverage/``.

The snapshots were audited against the original Python sincei 0.6.1 on
2026-09-21. The bins and the raw counts agree in every scenario, except for
these accepted differences:

- CPM: the original divides by the sum of its bedGraph lines after merging
  equal neighbouring bins, and skips a group with fewer than 10 counts; the port
  divides by the sum of all bin counts.
- The original ignores ``--minAlignedFraction`` and ``--ignoreForNormalization``.
- ``--duplicateFilter start_bc_umi`` keeps more reads in cluster B in the port;
  the port is correct.
- A group without signal: the original writes a bigwig without data, the port
  one zero-valued base. Both render to no bins.

Not covered: ``--motifFilter`` needs a 2bit genome that is not in ``testdata/``.

Run the tests::

    pytest tests/test_scBulkCoverage.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pybigtools
import pytest
from _cli_testing import (
    BAM1,
    BAM2,
    BAM_MERGED,
    DATA,
    Scenario,
    read_snapshot,
    run,
    tool_path,
)

if TYPE_CHECKING:
    from pathlib import Path

TOOL = "scBulkCoverage"
SL2 = DATA / "scCountReads_data"
FIXTURES = DATA / TOOL

BASE = [
    "-b",
    str(SL2 / "SL2-1.bam"),
    str(SL2 / "SL2-2.bam"),
    "--cellTag",
    "BC",
    "--numberOfProcessors",
    "1",
]
GROUP_INFO = ["--groupInfo", str(FIXTURES / "SL2_groups.tsv")]
# deeptools reads a region as chrom:start:end; the port accepts that form too.
REGION = ["--region", "chr1:23360000:23390000"]


@dataclass(frozen=True)
class BulkScenario(Scenario):
    """A :class:`Scenario` with the bin size its output is rendered at."""

    bin_size: int = 100
    whole_genome: bool = False


SCENARIOS: dict[str, BulkScenario] = {
    "baseline": BulkScenario(
        [], differs=False, reason="it is the reference the others are compared against"
    ),
    "normalize_none": BulkScenario(["--normalizeUsing", "None"]),
    "normalize_mean": BulkScenario(["--normalizeUsing", "Mean"]),
    "normalize_frequency": BulkScenario(["--normalizeUsing", "Frequency"]),
    "scale_factor": BulkScenario(["--normalizeUsing", "None", "--scaleFactor", "2.5"]),
    "ignore_for_normalization": BulkScenario(["--ignoreForNormalization", "chr1"]),
    "bin_size": BulkScenario(["--binSize", "1000"], bin_size=1000),
    "duplicate_filter": BulkScenario(["--duplicateFilter", "start_bc_umi"]),
    "extend_reads": BulkScenario(["--extendReads", "300"]),
    "center_reads": BulkScenario(["--extendReads", "300", "--centerReads"]),
    "min_mapping_quality": BulkScenario(["--minMappingQuality", "50"]),
    "sam_flag_exclude": BulkScenario(["--samFlagExclude", "16"]),
    "min_fragment_length": BulkScenario(["--minFragmentLength", "300"]),
    "max_fragment_length": BulkScenario(["--maxFragmentLength", "300"]),
    "filter_rna_strand": BulkScenario(["--filterRNAstrand", "forward"]),
    "min_aligned_fraction": BulkScenario(["--minAlignedFraction", "0.9"]),
    "gc_content": BulkScenario(["--GCcontentFilter", "0.30,0.50"]),
    "blacklist": BulkScenario(["--blacklist", str(SL2 / "test_regions.bed")]),
    "labels": BulkScenario(["--labels", "SL2-2", "SL2-1"]),
    "mnase": BulkScenario(["--mnase"]),
    "offset": BulkScenario(["--offset", "1"]),
    "whole_genome": BulkScenario(
        [],
        whole_genome=True,
        differs=False,
        reason="every read lies inside the baseline region",
    ),
}


def args_for(scenario: BulkScenario) -> list[str]:
    """The arguments of ``scenario``, without ``-o``."""
    region = [] if scenario.whole_genome else REGION
    return [*BASE, *GROUP_INFO, *region, *scenario.args]


def render(out_dir: Path, bin_size: int) -> str:
    """Project the bigwigs in ``out_dir`` to the text a reference can diff.

    The first line names every output file. Then one line per non-zero bin,
    ``cluster chrom start end value``, with runs of equal bins split back into
    bins of ``bin_size``: the original writes zero runs and merges equal
    neighbours, the port does neither, and neither choice carries meaning.
    """
    files = sorted(p.name for p in out_dir.iterdir())
    lines = ["files\t" + "\t".join(files)]
    for name in files:
        if not name.endswith(".bw"):
            continue
        cluster = name.removeprefix("cov_").removesuffix(".bw")
        with pybigtools.open(str(out_dir / name)) as bw:
            for chrom in sorted(bw.chroms()):
                for start, end, value in bw.records(chrom):
                    if value == 0:
                        continue
                    for bin_start in range(start, end, bin_size):
                        bin_end = min(bin_start + bin_size, end)
                        lines.append(
                            f"{cluster}\t{chrom}\t{bin_start}\t{bin_end}\t{value:.9g}"
                        )
    return "\n".join(lines) + "\n"


Coverage = dict[tuple[str, str, int, int], float]


def parse(text: str) -> tuple[str, Coverage]:
    """Split a rendering into its file line and its bins."""
    files, *rows = text.splitlines()
    bins: Coverage = {}
    for row in rows:
        cluster, chrom, start, end, value = row.split("\t")
        bins[cluster, chrom, int(start), int(end)] = float(value)
    return files, bins


def assert_same_coverage(got: str, want: str) -> None:
    """Same files and bins; values equal up to bigwig's float32 precision."""
    got_files, got_bins = parse(got)
    want_files, want_bins = parse(want)
    assert got_files == want_files
    missing = sorted(want_bins.keys() - got_bins.keys())
    extra = sorted(got_bins.keys() - want_bins.keys())
    assert not missing, f"{len(missing)} bins missing, first: {missing[:5]}"
    assert not extra, f"{len(extra)} bins extra, first: {extra[:5]}"
    keys = sorted(want_bins)
    np.testing.assert_allclose(
        [got_bins[k] for k in keys], [want_bins[k] for k in keys], rtol=1e-6
    )


def run_port(args: list[str], out_dir: Path, bin_size: int = 100) -> str:
    """Run the port into ``out_dir`` and return its rendered output."""
    out_dir.mkdir()
    proc = run(tool_path(TOOL), [*args, "-o", str(out_dir / "cov")])
    assert proc.returncode == 0, (
        f"{TOOL} {' '.join(args)} failed (exit {proc.returncode}):\n{proc.stderr}"
    )
    return render(out_dir, bin_size)


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_matches_its_snapshot(name: str, tmp_path: Path) -> None:
    scenario = SCENARIOS[name]
    got = run_port(args_for(scenario), tmp_path / "out", scenario.bin_size)
    assert_same_coverage(got, read_snapshot(TOOL, name))


def test_every_scenario_that_should_differ_does() -> None:
    """A scenario whose reference equals baseline asserts nothing about its flag."""
    baseline = read_snapshot(TOOL, "baseline")
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            continue
        assert read_snapshot(TOOL, name) != baseline, (
            f"scenario '{name}' ({' '.join(scenario.args)}) produces the baseline "
            f"output, so it cannot fail if the flag breaks. Pick a discriminating "
            f"value, or mark it differs=False with a reason."
        )


def test_scenarios_that_cannot_differ_say_why() -> None:
    for name, scenario in SCENARIOS.items():
        if not scenario.differs:
            assert scenario.reason, f"scenario '{name}' needs a recorded reason"


def test_a_group_without_reads_is_warned_about(tmp_path: Path) -> None:
    """The blacklist covers the locus, so neither group keeps a read."""
    proc = run(
        tool_path(TOOL),
        [*args_for(SCENARIOS["blacklist"]), "-o", str(tmp_path / "cov")],
    )

    assert proc.returncode == 0, proc.stderr
    for group in ("A", "B"):
        assert f'no reads were found for group "{group}"' in proc.stderr


# Properties that hold between two runs, without a snapshot.


def test_without_group_info_every_barcode_is_pooled_into_one_track(
    tmp_path: Path,
) -> None:
    """Every barcode in the SL2 BAMs is in one of the two clusters."""
    unnormalized = ["--normalizeUsing", "None"]
    _, grouped = parse(
        run_port([*BASE, *GROUP_INFO, *REGION, *unnormalized], tmp_path / "grouped")
    )
    files, pooled = parse(
        run_port([*BASE, *REGION, *unnormalized], tmp_path / "pooled")
    )

    assert files == "files\tcov.bw"
    summed: Coverage = {}
    for (_, chrom, start, end), value in grouped.items():
        key = ("cov", chrom, start, end)
        summed[key] = summed.get(key, 0.0) + value
    assert pooled.keys() == summed.keys()
    keys = sorted(summed)
    np.testing.assert_allclose(
        [pooled[k] for k in keys], [summed[k] for k in keys], rtol=1e-6
    )


@pytest.mark.parametrize("method", ["Mean", "Frequency"])
def test_without_group_info_cell_normalizations_are_rejected(
    method: str, tmp_path: Path
) -> None:
    proc = run(
        tool_path(TOOL),
        [*BASE, *REGION, "--normalizeUsing", method, "-o", str(tmp_path / "cov")],
    )

    assert proc.returncode != 0
    assert "need a group info file" in proc.stderr


def _bedgraph_bins(out_dir: Path) -> Coverage:
    bins: Coverage = {}
    for path in sorted(out_dir.glob("*.bedgraph")):
        cluster = path.stem.removeprefix("cov_")
        for line in path.read_text().splitlines():
            if line.startswith("track"):
                continue
            chrom, start, end, value = line.split("\t")
            if float(value) != 0:
                bins[cluster, chrom, int(start), int(end)] = float(value)
    return bins


def test_bedgraph_output_equals_the_bigwig(tmp_path: Path) -> None:
    args = args_for(SCENARIOS["baseline"])
    _, bigwig = parse(run_port(args, tmp_path / "bw"))

    out_dir = tmp_path / "bg"
    run_port([*args, "--outFileFormat", "bedgraph"], out_dir)
    bedgraph = _bedgraph_bins(out_dir)

    assert bedgraph.keys() == bigwig.keys()
    keys = sorted(bigwig)
    np.testing.assert_allclose(
        [bedgraph[k] for k in keys], [bigwig[k] for k in keys], rtol=1e-6
    )


def test_rpkm_is_cpm_per_kilobase_of_bin(tmp_path: Path) -> None:
    args = args_for(SCENARIOS["baseline"])
    _, cpm = parse(run_port(args, tmp_path / "cpm"))
    _, rpkm = parse(run_port([*args, "--normalizeUsing", "RPKM"], tmp_path / "rpkm"))

    assert rpkm.keys() == cpm.keys()
    keys = sorted(cpm)
    np.testing.assert_allclose(
        [rpkm[k] for k in keys], [cpm[k] * 1000 / 100 for k in keys], rtol=1e-6
    )


def test_a_merged_bam_grouped_by_tag_matches_its_separate_sources(
    tmp_path: Path,
) -> None:
    """One merged file split back apart by read group must equal its sources."""
    common = [
        "--groupInfo",
        str(DATA / "test_group_info.tsv"),
        "--cellTag",
        "BC",
        "--numberOfProcessors",
        "1",
        "--binSize",
        "10000",
    ]
    separate = run_port(["-b", BAM1, BAM2, *common], tmp_path / "separate", 10000)
    grouped = run_port(
        ["-b", BAM_MERGED, *common, "--groupTag", "RG"], tmp_path / "grouped", 10000
    )

    assert parse(separate)[1], "the sources produced no coverage at all"
    assert grouped == separate
