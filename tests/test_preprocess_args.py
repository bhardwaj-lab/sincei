#!/usr/bin/env python
"""Tests for ``preprocess_args``, the argv rewrite every command runs first.

Typer cannot express whitespace-separated multi-value options, so the CLI
rewrites ``--opt a b`` into ``--opt a --opt b`` before click ever sees it. The
rewrite has to tell an option apart from a value, and a negative number looks
like both.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from sincei.cli._parsers import preprocess_args

if TYPE_CHECKING:
    import pytest


def rewrite(argv: list[str], monkeypatch: pytest.MonkeyPatch) -> list[str]:
    monkeypatch.setattr(sys, "argv", ["sincei", "scBulkCoverage", *argv])
    preprocess_args()
    return sys.argv[2:]


def test_several_values_become_repeated_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert rewrite(["--offset", "5", "10"], monkeypatch) == [
        "--offset",
        "5",
        "--offset",
        "10",
    ]


def test_a_negative_value_stays_a_value(monkeypatch: pytest.MonkeyPatch) -> None:
    # `--offset 1 -1` is the documented way to say "the whole read". Treating
    # `-1` as an option left it stranded and click rejected it outright.
    assert rewrite(["--offset", "1", "-1"], monkeypatch) == [
        "--offset",
        "1",
        "--offset",
        "-1",
    ]


def test_a_lone_negative_value_survives(monkeypatch: pytest.MonkeyPatch) -> None:
    assert rewrite(["--offset", "-1"], monkeypatch) == ["--offset", "-1"]


def test_negative_values_on_both_sides(monkeypatch: pytest.MonkeyPatch) -> None:
    assert rewrite(["--offset", "-5", "-1"], monkeypatch) == [
        "--offset",
        "-5",
        "--offset",
        "-1",
    ]


def test_a_real_flag_after_a_value_still_starts_a_new_option(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Only a `-` followed by a digit is a value; every flag begins with a letter.
    assert rewrite(["--offset", "1", "-p", "4"], monkeypatch) == [
        "--offset",
        "1",
        "-p",
        "4",
    ]


def test_a_bare_extend_reads_gets_its_sentinel(monkeypatch: pytest.MonkeyPatch) -> None:
    # No value means "estimate from the data", which the backend reads as 0.
    assert rewrite(["--extendReads", "--centerReads"], monkeypatch) == [
        "--extendReads",
        "0",
        "--centerReads",
    ]


def test_a_boolean_flag_is_left_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    assert rewrite(["--mnase", "-p", "2"], monkeypatch) == ["--mnase", "-p", "2"]
