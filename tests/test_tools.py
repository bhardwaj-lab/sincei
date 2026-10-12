"""Every sincei console script starts and prints its help.

A command whose module cannot be imported exits before it does anything, and
the other suites do not run every command.
"""

from __future__ import annotations

from importlib.metadata import distribution

import pytest
from _cli_testing import run, tool_path

SCRIPTS = sorted(
    entry.name
    for entry in distribution("sincei").entry_points
    if entry.group == "console_scripts"
)


@pytest.mark.parametrize("script", SCRIPTS)
def test_help_exits_cleanly(script: str) -> None:
    proc = run(tool_path(script), ["-h"])
    assert proc.returncode == 0, (
        f"{script} -h failed (exit {proc.returncode}):\n{proc.stderr}"
    )
    assert "Usage" in proc.stdout
