"""Shared plumbing for the CLI snapshot suites.

Both `test_scFilterBarcodes.py` and `test_scFilterStats.py` drive an installed
console script, capture its TSV output and compare it against a committed
snapshot.  The pieces they share, locating the tool, running it and
normalizing its output, live here.

Snapshots are committed fixtures; nothing here writes them.  Regenerating one
means running the tool (or the reference implementation it came from) by hand
and replacing the file, so an output change is always a deliberate, reviewed
diff rather than a side effect of running the tests.

`Scenario` is the important one: it records not just the arguments to pass but
whether the resulting output is *expected to differ* from the baseline.  A
scenario whose snapshot equals the baseline asserts nothing about its flag, so
each suite has a guard test that fails when a scenario marked ``differs=True``
turns out to match.  Scenarios that legitimately cannot differ carry a
``reason`` instead.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "testdata"
BAM1 = str(DATA / "test_i1.bam")
BAM2 = str(DATA / "test_i2.bam")
BAM_MERGED = str(DATA / "test_i1_i2.bam")
BARCODES = str(DATA / "test_barcodes.txt")
BARCODES_1MIS = str(DATA / "test_barcodes_1mis.txt")
BED = str(DATA / "Chrna9_regions.bed")


@dataclass(frozen=True)
class Scenario:
    """One CLI invocation, and what its output is expected to do.

    ``differs`` is the contract: True means the flag must visibly change the
    output relative to ``baseline``.  When it cannot (because the flag only
    touches stderr, or the backend ignores it) set it False and say why in
    ``reason``, so the exception is a decision on the record rather than an
    oversight.
    """

    args: list[str] = field(default_factory=list)
    differs: bool = True
    reason: str = ""
    # Scenarios that set their own sampling flags opt out of the suite default.
    own_sampling: bool = False

    def __post_init__(self) -> None:
        if not self.differs and not self.reason:
            msg = "a scenario marked differs=False must record why"
            raise ValueError(msg)


def norm(text: str) -> str:
    """Normalize for comparison: strip trailing whitespace, drop blank tail.

    The original sincei writes empty tab-separated fields for all-zero rows
    where the port writes nothing, so trailing whitespace is not signal.
    """
    return "\n".join(line.rstrip() for line in text.splitlines()).rstrip("\n")


def tool_path(name: str) -> str:
    """Path to the installed port, or skip the test."""
    path = shutil.which(name)
    if path is None:
        pytest.skip(f"{name} not on PATH (run pytest inside the sincei_dev env)")
    return path


_ANSI = re.compile(r"\x1b\[[0-9;]*m")

# Typer prints its errors and tracebacks through rich, which picks colour and
# line width from the environment.  A runner that exports FORCE_COLOR makes
# rich colour each fragment separately, so substring assertions on the error
# text then fail in CI while passing locally.
# Pin the environment, and strip ANSI escapes.
_NO_COLOUR = {"NO_COLOR": "1", "TERM": "dumb", "COLUMNS": "200"}
_FORCED_COLOUR = ("FORCE_COLOR", "CLICOLOR_FORCE", "TTY_COMPATIBLE")


def run(tool: str, args: list[str]) -> subprocess.CompletedProcess[str]:
    """Run `tool` and hand back the completed process, success or not.

    Output comes back free of ANSI escapes, so tests can assert on the text a
    user reads rather than on how rich happened to decorate it.
    """
    env = {k: v for k, v in os.environ.items() if k not in _FORCED_COLOUR}
    env.update(_NO_COLOUR)
    proc = subprocess.run(
        [tool, *args], capture_output=True, text=True, check=False, env=env
    )
    return subprocess.CompletedProcess(
        proc.args,
        proc.returncode,
        _ANSI.sub("", proc.stdout),
        _ANSI.sub("", proc.stderr),
    )


def run_ok(tool: str, base: list[str], extra: list[str], out: str) -> str:
    """Run the tool, assert it succeeded, return the output file's text."""
    proc = run(tool, [*base, "-o", out, *extra])
    assert proc.returncode == 0, (
        f"{Path(tool).name} {' '.join(extra)} failed "
        f"(exit {proc.returncode}):\n{proc.stderr}"
    )
    return Path(out).read_text()


def snapshot_path(tool: str, name: str) -> Path:
    return DATA / tool / f"{name}.tsv"


def read_snapshot(tool: str, name: str) -> str:
    """Read a committed snapshot, failing loudly when it is missing.

    Deliberately a failure rather than a skip: a deleted fixture would
    otherwise silently reduce coverage without anything going red.
    """
    snap = snapshot_path(tool, name)
    if not snap.exists():
        pytest.fail(
            f"missing snapshot {snap}: regenerate it by running the tool with "
            f"this scenario's arguments and committing the output"
        )
    return snap.read_text()
