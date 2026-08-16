"""Parsers and shared plumbing between the argparse CLI and the Rust backend.

The CLI exposes user-friendly option shapes (choice strings, ``"low,high"``
strings, barcode/whitelist files) whereas the Rust functions take plain scalars
and lists.  This module centralizes those conversions and input checks.

This is the argparse frontend's own copy: it deliberately does not import from
``sincei.cli_typer``, because that package's ``__init__`` pulls in typer.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# CLI duplicate-filter choice -> Rust `dup_method` string.
DUPLICATE_FILTER_CHOICES = ["start_bc", "start_umi", "start_end_bc", "start_end_bc_umi"]

_DUP_METHOD_MAP: dict[str, str] = {
    "start_bc": "barcode_start",
    "start_umi": "barcode_umi_start",
    "start_end_bc": "barcode_start_end",
    "start_end_bc_umi": "barcode_umi_start_end",
}


def dup_method(value: str | None) -> str | None:
    """Map the CLI duplicate-filter choice to the backend's method string."""
    return _DUP_METHOD_MAP[value] if value is not None else None


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
                f"--labels count ({len(labels)}) does not match the number of "
                f"BAM files ({len(bam_files)})."
            )
            raise ValueError(msg)
        return list(labels)
    return [Path(b).stem for b in bam_files]


def warn_unsupported(**options: object) -> None:
    """Emit a warning for any CLI option that the backend does not yet honor."""
    for name, value in options.items():
        if value:
            logger.warning(
                "Option %r is not supported by this command and will be ignored.", name
            )


def log_parameters(**parameters: object) -> None:
    """Log each ``name: value`` parameter at INFO level (used under ``--verbose``)."""
    for name, value in parameters.items():
        logging.info("%s: %s", name, value)
