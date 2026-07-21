"""Typed wrappers for third-party decorators used across models.

These wrappers preserve runtime behavior while helping static analyzers
(mypy/Pylance) keep function signatures precise when decorator stubs are
insufficiently typed.
"""

from collections.abc import Callable
from typing import Any, Literal, TypeVar, cast

import pandera.pandas as pa
from pydantic import field_serializer, field_validator, model_validator

F = TypeVar("F", bound=Callable[..., Any])


def field_validator_typed(*field_names: str, **kwargs: Any) -> Callable[[F], F]:
    """Typed wrapper around pydantic ``field_validator``."""
    return cast(Callable[[F], F], field_validator(*field_names, **kwargs))


def field_serializer_typed(*field_names: str, **kwargs: Any) -> Callable[[F], F]:
    """Typed wrapper around pydantic ``field_serializer``."""
    return cast(Callable[[F], F], field_serializer(*field_names, **kwargs))


def check_types_typed(func: F) -> F:
    """Typed wrapper around pandera ``check_types``."""
    return pa.check_types(func)


def model_validator_typed(*, mode: Literal["wrap", "before", "after"] = "after") -> Callable[[F], F]:
    """Typed wrapper around pydantic ``model_validator``."""
    return cast(Callable[[F], F], model_validator(mode=mode))


def command_decorator_typed(decorator_factory: Callable[..., Any]) -> Callable[..., Callable[[F], F]]:
    """Typed wrapper for command decorator factories (e.g. ``Typer.command``)."""
    return cast(Callable[..., Callable[[F], F]], decorator_factory)
