from __future__ import annotations

import logging
import os
from collections.abc import Callable

import click
import typer

from sincei import _sincei as internal


def version_string() -> str:
    return internal.version()


def log_parameters(**parameters: object) -> None:
    for name, value in parameters.items():
        logging.info("%s: %s", name, value)


def normalize_region(ctx: click.Context, param: click.Parameter, string: str | None) -> str | None:
    if string is None:
        return None

    # Remove whitespaces using split,join trick
    region = "".join(string.split())
    if region == "":
        return None

    # Remove undesired characters that may be present and replace "-" by ":"
    region = region.translate({ord(character): None for character in ",;|!{}()"}).replace("-", ":")

    if len(region) == 0:
        raise click.BadParameter(f"{string} is not a valid region", param=param, ctx=ctx)

    return region


def normalize_processors(ctx: click.Context, param: click.Parameter, value: str) -> int:
    # Use os.sched_getaffinity if available (Linux)
    # Fallback to os.cpu_count()
    try:
        available_processors = len(os.sched_getaffinity(0))
    except AttributeError:
        available_processors = os.cpu_count() or 1

    if value == "max/2":
        return int(available_processors * 0.5)
    if value == "max":
        return available_processors

    try:
        number_of_processors = int(value)
    except ValueError as exc:
        raise click.BadParameter(f"{value} is not a valid number of processors", param=param, ctx=ctx) from exc

    if number_of_processors > available_processors:
        number_of_processors = available_processors

    return number_of_processors


def smart_label(label: str) -> str:
    basename = os.path.splitext(os.path.basename(label))[0]
    if basename == "":
        basename = os.path.basename(label)
    return basename


def smart_labels(labels: list[str]) -> list[str]:
    inferred = [smart_label(label) for label in labels]
    if len(inferred) != len(set(inferred)):
        print(
            "Labels inferred from file names are not unique. "
            "Please be aware that in case of overlapping barcodes the counts will be merged."
        )
    return inferred


def version_option(prog_name: str) -> Callable[[Callable[..., object]], Callable[..., object]]:
    def decorator(function: Callable[..., object]) -> Callable[..., object]:
        def callback(ctx: click.Context, param: click.Parameter, value: bool) -> bool:
            if value:
                typer.echo(f"{prog_name} {version_string()}")
                raise typer.Exit()
            return value

        return click.option(
            "-V",
            "--version",
            is_flag=True,
            is_eager=True,
            expose_value=False,
            help="Print the program version and exit.",
            callback=callback,
        )(function)

    return decorator
