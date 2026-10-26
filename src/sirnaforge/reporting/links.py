"""Deterministic, offline source resolvers for explicitly named biological identifiers."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any
from urllib.parse import quote, urlsplit


def safe_source_url(value: str) -> bool:
    """Allow only HTTPS source pages on Ensembl or MirGeneDB, without credentials or active content."""
    try:
        parsed = urlsplit(value)
        host = parsed.hostname or ""
        valid_host = host in {"www.ensembl.org", "ensembl.org", "mirgenedb.org", "www.mirgenedb.org"} or bool(
            re.fullmatch(r"(?:e\d+|[a-z]{3}\d{4}\.archive)\.ensembl\.org", host)
        )
        return (
            parsed.scheme == "https"
            and valid_host
            and parsed.username is None
            and parsed.password is None
            and parsed.port is None
            and not any(c in value for c in "\\<>\"'\r\n\t")
            and not re.search(r"%(?:0[ad]|22|27|3[ce]|5c)", value, re.IGNORECASE)
        )
    except ValueError:
        return False


def source_link(
    identifier: str, source: str, *, release: str | None = None, resolver: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Resolve a supported identifier locally, leaving unsupported identities explicitly unresolved."""
    result: dict[str, Any] = {"identifier": identifier, "source": source, "release": release, "url": None}
    if source == "ensembl" and re.fullmatch(r"ENS[A-Z]*[GT]\d+(?:\.\d+)?", identifier):
        template = (
            "https://e{release}.ensembl.org/id/{identifier}"
            if release and release.isdigit()
            else "https://www.ensembl.org/id/{identifier}"
        )
        values = {"identifier": quote(identifier, safe=""), "release": quote(release or "", safe="")}
    elif source == "mirgenedb" and re.fullmatch(r"[A-Z][a-z]{2}-Mir-[A-Za-z0-9-]+(?:_[35]p\*?)?", identifier):
        species, name = identifier.split("-", 1)
        name = re.sub(r"_[35]p\*?$", "", name)
        template = "https://mirgenedb.org/show/{species}/{name}"
        values = {
            "species": species.lower(),
            "name": quote(name, safe=""),
            "identifier": quote(identifier, safe=""),
            "release": quote(release or "", safe=""),
        }
    else:
        result["reason"] = "identifier_not_supported"
        return result
    resolver_data = resolver or {}
    declared = (
        resolver_data.get("release_template")
        if source == "ensembl" and release and release.isdigit()
        else resolver_data.get("template")
    )
    if declared is not None:
        if not isinstance(declared, str):
            result["reason"] = "invalid_resolver_template"
            return result
        template = declared
    try:
        url = template.format_map(values)
    except (KeyError, ValueError, AttributeError):
        result["reason"] = "invalid_resolver_template"
        return result
    if not safe_source_url(url):
        result["reason"] = "unsafe_resolver_url"
        return result
    result.update(url=url, template=template, reason=None)
    return result


def canonical_identity_link(
    identifier: str,
    canonical_source: str | None,
    canonical_release: str | None,
    *,
    resolver: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Link an Ensembl ID without treating another authority's release as Ensembl's."""
    ensembl_release = None if canonical_source and canonical_source.startswith("GENCODE") else canonical_release
    return source_link(identifier, "ensembl", release=ensembl_release, resolver=resolver)
