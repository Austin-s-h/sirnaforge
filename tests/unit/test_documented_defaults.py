"""Documented defaults must be the defaults the code applies (#99 acceptance).

Two documented numbers had drifted from the code by `8dce4ae`: the off-target cap (documented 3,
enforced 15) and the design weights (documented 0.40/0.35, applied 0.35/0.40). Restating a number in
prose is what allows that, so these tests read the number back out of the CLI help and the
documentation tables and compare it against the model defaults. They are the mechanism that makes
the fix stick, not a record that it was made once.
"""

import re
from pathlib import Path
from typing import get_args

import pytest
import typer

from sirnaforge.cli import app
from sirnaforge.config.run_policy import default_for, mirna_preset_default_for
from sirnaforge.models.sirna import (
    DesignParameters,
    DesignWeights,
    FilterCriteria,
    MiRNADesignConfig,
    OffTargetFilterCriteria,
    PostScreenMiRNAWeights,
    PostScreenSiRNAWeights,
    SiRNACandidate,
    TargetAccessibilityConfig,
)

REPO_ROOT = Path(__file__).resolve().parents[2]

# Options whose help text quotes a default, and the policy setting that default comes from.
HELP_DEFAULTS: tuple[tuple[str, str, str], ...] = (
    ("workflow", "gc_min", "gc_min"),
    ("workflow", "gc_max", "gc_max"),
    ("workflow", "sirna_length", "sirna_length"),
    ("workflow", "max_off_targets", "max_off_target_count"),
    ("workflow", "min_asymmetry", "min_asymmetry_score"),
    ("workflow", "max_paired_fraction", "max_paired_fraction"),
    ("workflow", "plfold_window", "plfold_window"),
    ("workflow", "plfold_max_bp_span", "plfold_max_bp_span"),
    ("workflow", "accessibility_log_floor", "accessibility_log_floor"),
    ("design", "gc_min", "gc_min"),
    ("design", "gc_max", "gc_max"),
    ("design", "length", "sirna_length"),
    ("design", "max_poly_runs", "max_poly_runs"),
    ("design", "min_asymmetry", "min_asymmetry_score"),
    ("design", "max_paired_fraction", "max_paired_fraction"),
    ("design", "min_empirical", "min_empirical_score"),
    ("design", "plfold_window", "plfold_window"),
    ("design", "plfold_max_bp_span", "plfold_max_bp_span"),
    ("design", "accessibility_log_floor", "accessibility_log_floor"),
)

# Model classes the documentation quotes field defaults from, by the class name in the code block.
DOCUMENTED_MODELS = {
    "FilterCriteria": FilterCriteria,
    "OffTargetFilterCriteria": OffTargetFilterCriteria,
    "TargetAccessibilityConfig": TargetAccessibilityConfig,
    "DesignParameters": DesignParameters,
    "MiRNADesignConfig": MiRNADesignConfig,
    "DesignWeights": DesignWeights,
    "PostScreenSiRNAWeights": PostScreenSiRNAWeights,
    "PostScreenMiRNAWeights": PostScreenMiRNAWeights,
}

CLASS_LINE = re.compile(r"^class ([A-Za-z]+)\(")
FIELD_LINE = re.compile(r"^ {4}([a-z_0-9]+):\s*([A-Za-z|\[\]0-9 ]+?)\s*=\s*([^#\n]+?)\s*(?:#.*)?$")
DEFAULT_IN_HELP = re.compile(r"default:\s*(-?[0-9]+(?:\.[0-9]+)?)")


def _option_help(command_name: str, parameter: str) -> str:
    """The help string of one option, taken from the Click command rather than rendered output."""
    command = typer.main.get_command(app)
    subcommand = command.commands[command_name]  # type: ignore[attr-defined]
    for param in subcommand.params:
        if param.name == parameter:
            return str(param.help or "")
    raise AssertionError(f"{command_name} has no parameter named {parameter!r}")


def _documented_fields(path: Path):
    """Yield ``(line_number, class_name, field, annotation, literal)`` for documented field defaults."""
    in_block = False
    current: str | None = None
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if line.startswith("```"):
            in_block = line.startswith("```python")
            current = None
            continue
        if not in_block:
            continue
        class_match = CLASS_LINE.match(line)
        if class_match:
            current = class_match.group(1)
            continue
        field_match = FIELD_LINE.match(line)
        if field_match and current in DOCUMENTED_MODELS:
            yield number, current, field_match.group(1), field_match.group(2).strip(), field_match.group(3).strip()


def _as_value(literal: str):
    """Parse a documented literal into the value it claims, or ``None`` if it is not one."""
    text = literal.strip().rstrip(",")
    if text in {"None", "True", "False"}:
        return {"None": None, "True": True, "False": False}[text]
    if text.startswith(('"', "'")):
        return text.strip("\"'")
    try:
        return float(text) if "." in text else int(text)
    except ValueError:
        return None


@pytest.mark.unit
def test_every_default_quoted_in_cli_help_matches_the_resolved_default():
    """A help string that names a number must name the one the resolver applies.

    `--max-off-targets` said "default: 3" while `OffTargetFilterCriteria` enforced 15, so a user
    setting 10 to tighten the gate was loosening it.
    """
    for command_name, parameter, setting in HELP_DEFAULTS:
        help_text = _option_help(command_name, parameter)
        quoted = DEFAULT_IN_HELP.search(help_text)
        assert quoted, f"{command_name} --{parameter} quotes no default; it should name one"
        assert float(quoted.group(1)) == float(default_for(setting)), (
            f"{command_name} --{parameter} help says {quoted.group(1)}, the resolver applies {default_for(setting)}"
        )


@pytest.mark.unit
def test_the_options_the_mirna_preset_moves_name_both_numbers():
    """A miRNA user reading only the siRNA default reads a number that never applies to them."""
    for command_name in ("workflow", "design"):
        for parameter, setting in (("gc_min", "gc_min"), ("gc_max", "gc_max"), ("overhang", "default_overhang")):
            help_text = _option_help(command_name, parameter)
            preset = mirna_preset_default_for(setting)
            if preset == default_for(setting):
                continue
            assert "mirna" in help_text, f"{command_name} --{parameter} does not mention the miRNA default"
            assert str(preset) in help_text


@pytest.mark.unit
def test_field_defaults_quoted_in_the_documentation_match_the_models():
    """The parameter tables are checked against the models rather than maintained beside them.

    Only ``models_and_scoring.md`` quotes field defaults as code; ``scoring.md`` states the weight
    vectors as a table, which the weight test below covers.
    """
    document = "docs/models_and_scoring.md"
    path = REPO_ROOT / document
    checked = 0
    for number, class_name, field, _annotation, literal in _documented_fields(path):
        model = DOCUMENTED_MODELS[class_name]
        if field not in model.model_fields:
            continue
        documented = _as_value(literal)
        if documented is None and literal.strip().rstrip(",") != "None":
            continue
        actual = model.model_fields[field].get_default(call_default_factory=False)
        checked += 1
        assert documented == pytest.approx(actual) if isinstance(actual, float) else documented == actual, (
            f"{document}:{number} documents {class_name}.{field} = {literal}, the model default is {actual}"
        )
    assert checked, f"{document} quoted no model field defaults; the parser has stopped matching"


@pytest.mark.unit
def test_documented_field_annotations_agree_with_the_models_on_nullability():
    """A default checker that reads only the literal lets the *type* drift, and here it matters.

    ``--filter-action <id>=off`` works by clearing a threshold to ``None``, so whether a threshold is
    nullable is a documented capability, not a detail. The docs said ``max_off_target_count: int = 15``
    where the model is ``int | None``.
    """
    path = REPO_ROOT / "docs/models_and_scoring.md"
    checked = 0
    for number, class_name, field, annotation, _literal in _documented_fields(path):
        model = DOCUMENTED_MODELS[class_name]
        if field not in model.model_fields:
            continue
        documented_nullable = "None" in annotation
        actual_nullable = type(None) in get_args(model.model_fields[field].annotation)
        checked += 1
        assert documented_nullable == actual_nullable, (
            f"docs/models_and_scoring.md:{number} documents {class_name}.{field} as {annotation!r}; "
            f"the model field is {'nullable' if actual_nullable else 'not nullable'}"
        )
    assert checked, "no documented annotations were compared; the parser has stopped matching"


@pytest.mark.unit
def test_the_documented_off_target_cap_is_the_enforced_one():
    """Both halves of the `8dce4ae` cap drift: the field description and the documentation."""
    cap = OffTargetFilterCriteria.model_fields["max_off_target_count"].get_default(call_default_factory=False)
    assert cap == 15

    count_description = str(SiRNACandidate.model_fields["off_target_count"].description)
    assert "≤3" not in count_description, "the reported column must not advertise a ceiling nothing applies"
    assert "max_off_target_count" in count_description

    for document in ("docs/models_and_scoring.md",):
        text = (REPO_ROOT / document).read_text()
        assert "max_off_target_count: int | None = 15" in text


@pytest.mark.unit
def test_the_documented_design_weights_are_the_applied_ones():
    """`design_v4` applies asymmetry 0.40 / accessibility 0.35; the docs said the reverse."""
    weights = DesignWeights()
    assert (weights.asymmetry, weights.target_accessibility) == (0.40, 0.35)

    for document in ("docs/models_and_scoring.md", "docs/scoring.md"):
        text = (REPO_ROOT / document).read_text()
        assert "target_accessibility  0.40" not in text
        assert "target_accessibility: float = 0.40" not in text
