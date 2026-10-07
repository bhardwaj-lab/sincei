"""Render a typer command line tool in the style sphinx-argparse used.

``.. typer-cli:: sincei.cli.scFilterBarcodes`` writes the command's description,
a usage line, and one section per help panel. Each option is a definition list
entry: its names, its help, and its default or its possible choices. A group
such as ``scCountReads`` gets one section per subcommand.

sphinx-argparse cannot read typer, and sphinx-click does not accept the command
classes typer builds on its own copy of click, hence this directive.
"""

from __future__ import annotations

import enum
import importlib
import re
from typing import TYPE_CHECKING, Any

import typer
from docutils import nodes
from docutils.parsers.rst import directives, roles
from docutils.statemachine import StringList
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles

if TYPE_CHECKING:
    from sphinx.application import Sphinx


# Rich markup such as ``[bold yellow]gzip[/bold yellow]``, which the terminal
# highlights. The help texts use it for values, so it becomes inline code.
_MARKUP = re.compile(r"\[(?P<style>[a-z][a-z ]*)\](?P<text>.*?)\[/(?P=style)\]", re.S)
# A paragraph that opens with a highlighted value and a colon describes that
# value, as in ``none: no compression``. Such paragraphs become list items.
_VALUE_ITEM = re.compile(r"``[^`]+``\s*:")


def _plain(help_text: str | None) -> str:
    """Help text as reST, with rich markup turned into inline code."""
    text = _MARKUP.sub(lambda m: f"``{m['text']}``", help_text or "")
    paragraphs = []
    for paragraph in text.split("\n\n"):
        if _VALUE_ITEM.match(paragraph.strip()):
            item = paragraph.strip().splitlines()
            paragraph = "\n".join(["* " + item[0], *("  " + line for line in item[1:])])
        paragraphs.append(paragraph)
    text = "\n\n".join(paragraphs)
    # A bullet list needs a blank line before it in reST, not in the terminal.
    lines: list[str] = []
    for line in text.splitlines():
        if (
            line.lstrip().startswith(("* ", "- "))
            and lines
            and lines[-1].strip()
            and (not lines[-1].lstrip().startswith(("* ", "- ")))
        ):
            lines.append("")
        lines.append(line)
    return "\n".join(lines)


def _value(value: Any) -> str:
    return str(value.value if isinstance(value, enum.Enum) else value)


def _names(param: Any) -> str:
    names = [*param.opts, *param.secondary_opts]
    return ", ".join(sorted(names, key=lambda name: not name.startswith("--")))


def _badge(role: str, text: str) -> str:
    escaped = text.replace("\\", "\\\\").replace("`", "\\`")
    return f":{role}:`{escaped}`"


def _badges(param: Any) -> list[str]:
    """The marks that close an option's help: required, or its default."""
    if param.required:
        return [_badge("cli-required", "required")]
    # A help text that states its own default (a computed or conditional one)
    # is left to say it.
    if param.show_default is False or "Default:" in (param.help or ""):
        return []
    if isinstance(param.show_default, str):
        return [_badge("cli-default", f"default: {param.show_default}")]
    if param.is_flag:
        return [_badge("cli-default", "default: True")] if param.default else []
    if param.default in (None, (), []):
        return []
    return [_badge("cli-default", f"default: {_value(param.default)}")]


def _ends_in_list(lines: list[str]) -> bool:
    """Whether the last paragraph of ``lines`` is a list item."""
    paragraph: list[str] = []
    for line in reversed(lines):
        if not line.strip():
            break
        paragraph.append(line)
    return bool(paragraph) and paragraph[-1].lstrip().startswith(("* ", "- "))


def _choices(param: Any) -> str | None:
    choices = getattr(param.type, "choices", None)
    # A help text that already names every choice says it better.
    named = _plain(param.help)
    if not choices or all(f"``{choice}``" in named for choice in choices):
        return None
    return "Possible choices: " + ", ".join(f"``{c}``" for c in choices)


def _usage(prog: str, command: Any) -> str:
    required = [
        f"{param.opts[-1]} {param.metavar or param.name.upper()}"
        for param in command.params
        if param.param_type_name == "option" and param.required
    ]
    return " ".join(["usage:", prog, *required, "[OPTIONS]"])


def _command_lines(prog: str, command: Any, underline: str) -> list[str]:
    lines = [*_plain(command.help).splitlines(), "", ".. code-block:: text", ""]
    lines += ["   " + _usage(prog, command), ""]

    panels: dict[str, list[Any]] = {}
    for param in command.params:
        if param.param_type_name != "option" or getattr(param, "hidden", False):
            continue
        panel = (getattr(param, "rich_help_panel", None) or "Options").strip()
        panels.setdefault(panel, []).append(param)

    for panel, params in panels.items():
        lines += [panel, underline * len(panel), ""]
        for param in params:
            lines.append(f"``{_names(param)}``")
            body = [*_plain(param.help).splitlines()]
            badges = " ".join(_badges(param))
            # The marks close the help text's last line, unless that line ends
            # a list, where they would read as part of its last item.
            if badges and body and not _ends_in_list(body):
                body[-1] = f"{body[-1].rstrip()} {badges}"
            elif badges:
                body += ["", badges]
            choices = _choices(param)
            if choices:
                body += ["", choices]
            lines += ["   " + line if line.strip() else "" for line in body]
            lines.append("")
    return lines


class TyperCli(SphinxDirective):
    """``.. typer-cli:: module.path`` with ``:prog:`` naming the program."""

    required_arguments = 1
    option_spec = {"prog": directives.unchanged_required}

    def run(self) -> list[nodes.Node]:
        module = importlib.import_module(self.arguments[0])
        command = typer.main.get_command(module.app)
        prog = self.options["prog"]

        subcommands = getattr(command, "commands", None) or {}
        if subcommands:
            lines = [*_plain(command.help).splitlines(), ""]
            for name, subcommand in subcommands.items():
                title = f"{prog} {name}"
                lines += [title, "-" * len(title), ""]
                lines += _command_lines(title, subcommand, "~")
        else:
            lines = _command_lines(prog, command, "-")

        container = nodes.section()
        container.document = self.state.document
        nested_parse_with_titles(
            self.state, StringList(lines, source=self.arguments[0]), container
        )
        return container.children


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_directive("typer-cli", TyperCli)
    # Inline marks, styled in content/_static/custom.css.
    for role in ("cli-required", "cli-default"):
        app.add_role(
            role, roles.CustomRole(role, roles.generic_custom_role, {"class": [role]})
        )
    return {"parallel_read_safe": True, "parallel_write_safe": True}
