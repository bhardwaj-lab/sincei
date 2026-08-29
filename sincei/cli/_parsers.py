"""Parsers and shared plumbing for the typer CLI and the Rust (`_sincei`) backend.

The CLI exposes user-friendly option shapes (choice enums, ``"low,high"``
strings, barcode/whitelist files) whereas the Rust functions take plain
scalars and lists.  This module centralizes those conversions and input
checks, plus the small amount of app plumbing shared by every command
(``configure_logging``, ``preprocess_args`` and ``log_parameters``).
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import typer

if TYPE_CHECKING:
    from collections.abc import Sequence

    import pandas as pd
    from anndata import AnnData

    from ._common_args import DuplicateFilter

logger = logging.getLogger(__name__)

_DUP_METHOD_MAP: dict[str, str] = {
    "start_bc": "barcode_start",
    "start_bc_umi": "barcode_umi_start",
    "start_end_bc": "barcode_start_end",
    "start_end_bc_umi": "barcode_umi_start_end",
}


def dup_method(value: DuplicateFilter | None) -> str | None:
    """Map the CLI duplicate-filter choice to the backend's method string."""
    return _DUP_METHOD_MAP[value.value] if value is not None else None


def umi_tag_if_used(
    tag: str | None, duplicate_filter: DuplicateFilter | None
) -> str | None:
    """The UMI tag to ask the backend for, or ``None`` when nothing reads it.

    Only the UMI-aware duplicate filters look at a record's UMI. Requesting the
    tag anyway would make every read pay a scan of its auxiliary-tag block for a
    value that is then discarded, on every command whether or not it
    deduplicates at all.
    """
    if tag is None or duplicate_filter is None:
        return None
    return tag if "umi" in duplicate_filter.value else None


def parse_motif_filter(motifs: list[str] | None) -> list[tuple[str, str]] | None:
    """Parse ``["A,TA", ...]`` into ``[("A", "TA"), ...]`` (read, reference)."""
    if not motifs:
        return None
    parsed: list[tuple[str, str]] = []
    for entry in motifs:
        parts = entry.split(",")
        if len(parts) != 2:
            msg = f"--motifFilter expects 'read_motif,ref_motif' (got {entry!r})"
            raise ValueError(msg)
        parsed.append((parts[0].strip(), parts[1].strip()))
    return parsed


def parse_gc_content(value: str | None) -> tuple[float | None, float | None]:
    """Parse ``"<low>,<high>"`` into ``(min_gc, max_gc)``; ``(None, None)`` if unset."""
    if not value:
        return None, None
    parts = value.split(",")
    if len(parts) != 2:
        msg = f"--GCcontentFilter expects '<low>,<high>' (got {value!r})"
        raise ValueError(msg)
    return float(parts[0]), float(parts[1])


def optional_length(value: int | None) -> int | None:
    """CLI fragment-length options use 0 to mean "no limit"; map that to ``None``."""
    return value if value and value > 0 else None


def first_blacklist(blacklist: list[str] | None) -> str | None:
    """The backend accepts a single blacklist file; use the first one provided."""
    if not blacklist:
        return None
    if len(blacklist) > 1:
        logger.warning(
            "Multiple blacklist files provided; only the first (%s) is used.",
            blacklist[0],
        )
    return blacklist[0]


def read_barcodes(path: str) -> list[str]:
    """Read a single-column barcode/whitelist file into a list of barcodes."""
    with Path(path).open() as handle:
        return [line.strip() for line in handle if line.strip()]


def resolve_labels(
    bam_files: list[str],
    labels: list[str] | None,
    use_smart_labels: bool,
) -> list[str]:
    """Resolve the per-BAM sample labels used to match the group-info file.

    Priority: explicit ``--labels`` (must match the BAM count), then
    ``--smartLabels`` / default, both of which use the file stem.
    """
    if labels:
        if len(labels) != len(bam_files):
            msg = (
                f"--labels count ({len(labels)}) does not match the number of BAM "
                f"files ({len(bam_files)}). Give one label per BAM, or omit "
                f"--labels to use the file names."
            )
            raise typer.BadParameter(msg)
        return list(labels)
    return [Path(b).stem for b in bam_files]


def require_single_bam_for_group_tag(
    bam_files: list[str],
    group_tag: str | None,
) -> None:
    """Reject more than one input BAM when ``--groupTag`` is set.

    A merged BAM carries each read's sample of origin in a tag, which is the
    only reason the option exists; several files would mean two competing
    notions of what a sample is.  The tools that accept ``--groupTag`` all need
    this, so it lives here rather than being restated in each.
    """
    if group_tag is not None and len(bam_files) > 1:
        msg = (
            f"--groupTag expects a single merged BAM, but {len(bam_files)} were "
            f"given. Merge them first, e.g. `samtools merge -r`."
        )
        raise typer.BadParameter(msg)


def warn_labels_ignored_under_group_tag(
    labels: list[str] | None,
    use_smart_labels: bool,
) -> None:
    """Warn that the label options have nothing to name under ``--groupTag``."""
    given = [
        flag
        for flag, passed in (("--labels", labels), ("--smartLabels", use_smart_labels))
        if passed
    ]
    if given:
        logger.warning(
            "%s ignored: with --groupTag the sample names come from the BAM's "
            "@RG read groups, not from the input file names.",
            " and ".join(given),
        )


def warn_unsupported(**options: object) -> None:
    """Emit a warning for any CLI option that the backend does not yet honor."""
    for name, value in options.items():
        if value:
            logger.warning(
                "Option %r is not supported by this command and will be ignored.", name
            )


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
# Columns every sincei-generated .h5ad is expected to carry, and their type.

_REQUIRED_OBS: dict[str, object] = {
    "sample": "category",
    "barcode": "category",
}
_REQUIRED_VAR: dict[str, object] = {
    "chrom": "category",
    "name": "category",
    "start": int,
    "end": int,
}


def _coerce_column(
    frame: pd.DataFrame, column: str, dtype: object, label: str
) -> str | None:
    """Coerce one column in place, returning a problem description or None."""
    if column not in frame.columns:
        return f"  .{label} is missing the column '{column}'"
    try:
        if frame[column].dtype != dtype:
            frame[column] = frame[column].astype(dtype)
    except (TypeError, ValueError) as exc:
        return f"  .{label}['{column}'] cannot be read as {dtype}: {exc}"
    return None


def _anndata_problems(adata: AnnData, filename: str | None = None) -> str | None:
    """Describe why `adata` is not valid sincei input, or None if it is."""
    problems: list[str | None] = []
    if adata.n_obs == 0:
        problems.append("  .obs has 0 observations (cells)")
    if adata.n_vars == 0:
        problems.append("  .var has 0 variables (features)")
    for column, dtype in _REQUIRED_OBS.items():
        problems.append(_coerce_column(adata.obs, column, dtype, "obs"))  # ty: ignore[invalid-argument-type]
    for column, dtype in _REQUIRED_VAR.items():
        problems.append(_coerce_column(adata.var, column, dtype, "var"))  # ty: ignore[invalid-argument-type]

    found = [problem for problem in problems if problem is not None]
    if not found:
        return None

    where = f"'{filename}'" if filename else "the input .h5ad file"
    return "\n".join([f"{where} is not a valid sincei-formatted input file:", *found])


def validate_anndata(adata: AnnData, filename: str | None = None) -> AnnData:
    """Check that an AnnData carries the fields the sincei tools require.

    Returns the same object with `.var["start"]` / `.var["end"]` cast to int.

    Exits with status 1 after describing every problem found, rather than
    raising, so a user gets one readable message instead of a traceback.
    """
    problems = _anndata_problems(adata, filename)
    if problems is not None:
        sys.stderr.write(problems + "\n")
        raise typer.Exit(code=1)
    return adata


def validate_anndata_list(
    adatas: Sequence[AnnData], filenames: Sequence[str] | None = None
) -> list[AnnData]:
    """`validate_anndata` for several inputs, reporting every bad file at once."""
    names = list(filenames) if filenames is not None else [None] * len(adatas)
    problems = [
        problem
        for adata, name in zip(adatas, names, strict=False)
        if (problem := _anndata_problems(adata, name)) is not None
    ]
    if problems:
        sys.stderr.write("\n".join(problems) + "\n")
        raise typer.Exit(code=1)
    return list(adatas)


_REGION = re.compile(r"([^:]+):(\d+)(?:-|:)(\d+)")


def parse_region(region: str) -> tuple[str, int | None, int | None]:
    """Split `CHROM`, `CHROM:START-END` or `CHROM:START:END` into its parts.

    A chromosome on its own yields `(chrom, None, None)`, meaning the whole
    chromosome. Complements `normalize_region`, which tidies the string a user
    typed; this turns the tidied string into coordinates.
    """
    text = region.strip()
    if ":" not in text:
        return text, None, None

    match = _REGION.fullmatch(text)
    if match is None:
        msg = (
            f"Invalid region {region!r}. Expected 'CHROM', 'CHROM:START-END' "
            f"or 'CHROM:START:END'."
        )
        raise typer.BadParameter(msg)

    chrom, start_text, end_text = match.groups()
    start, end = int(start_text), int(end_text)
    if start > end:
        msg = f"Invalid region {region!r}. START must be <= END."
        raise typer.BadParameter(msg)
    return chrom, start, end


# ---------------------------------------------------------------------------
# App plumbing shared by every command
# ---------------------------------------------------------------------------
def log_parameters(**parameters: object) -> None:
    """Log each ``name: value`` parameter at INFO level (used under ``--verbose``)."""
    for name, value in parameters.items():
        logging.info("%s: %s", name, value)


def configure_logging() -> None:
    """Set up logging for a CLI run.

    The library itself configures nothing, so that importing sincei leaves the
    root logger of the importing application alone. Every command calls this
    from its ``cli()`` entry point instead. ``basicConfig`` does nothing if the
    root logger already has a handler, so an application that wraps a command
    keeps its own setup.
    """
    logging.basicConfig(
        stream=sys.stderr,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )


_NEGATIVE_NUMBER = re.compile(r"^-\d|^-\.\d")


def _is_option(token: str) -> bool:
    """Whether an argv token names an option rather than being a value.

    A leading ``-`` is not enough: ``--offset 5 -1`` ends in a negative
    *number*, and treating it as an option leaves it stranded for click to
    reject with "No such option: -1". Every flag here begins with a letter, so
    a ``-`` followed by a digit is always a value. `argparse` applies the same
    rule, which is why the reference implementation accepts the form.
    """
    return token.startswith("-") and not _NEGATIVE_NUMBER.match(token)


def preprocess_args() -> None:
    """Expand whitespace-separated multi-value options into repeated options.

    Typer does not support whitespace-separated multi-value options, so we
    preprocess ``sys.argv`` so that::

        app some_command --filters f1 f2 f3 --envs e1 e2

    becomes::

        app some_command --filters f1 --filters f2 --filters f3 --envs e1 --envs e2

    //!\\ DOWNSIDE: options must always come after positional arguments. //!\\
    //!\\ DOWNSIDE: options must either take one value or repeated values. //!\\
    This is fine for sincei since all commands only use options.
    """

    # `--extendReads` may be given with no value (meaning "estimate from data").
    _EXTEND_READS_FLAGS = {"-e", "--extendReads"}

    # Collect the leading command path (everything up to the first option).
    final_cmd: list[str] = []
    for arg in sys.argv:
        if _is_option(arg):
            break
        final_cmd.append(arg)

    # Expand each option so every value gets its own flag.
    for idx, arg in enumerate(sys.argv):
        if _is_option(arg):
            opt_values: list[str] = []
            for value in sys.argv[idx + 1 :]:
                if _is_option(value):
                    break
                opt_values.append(value)

            if len(opt_values) >= 1:
                for opt_value in opt_values:
                    final_cmd.extend([arg, opt_value])
            elif arg in _EXTEND_READS_FLAGS:
                # A bare --extendReads (no value) means "estimate the fragment
                # length from the data". Typer cannot express an option with an
                # optional value, so we pass a "0" that the Rust backend
                # interprets as "estimate".
                final_cmd.extend([arg, "0"])
            else:
                final_cmd.append(arg)

    sys.argv = final_cmd
