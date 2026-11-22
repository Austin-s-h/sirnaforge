"""Resource boundary assertions shared by report fixtures."""

from __future__ import annotations

import re
from html.parser import HTMLParser

from sirnaforge.reporting.links import safe_source_url


class _Resources(HTMLParser):
    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        """Allow source anchors, while refusing remote assets and other active resources."""
        attributes = dict(attrs)
        for key in ("src", "href", "action", "poster", "data", "srcset"):
            value = attributes.get(key)
            if value and not value.startswith("#"):
                assert tag == "a" and key == "href" and safe_source_url(value), (tag, key, value)
                assert attributes.get("target") == "_blank"
                assert set((attributes.get("rel") or "").split()) >= {"noopener", "noreferrer"}
        assert tag not in {"iframe", "object", "embed", "base"}


def assert_offline_report(html: str) -> None:
    """Allow deterministic source links, never external scripts, assets, fetches or storage."""
    _Resources().feed(html)
    for url in re.findall(r"https?://[^\s\"'<>]+", html):
        # Resolver templates are recorded as data as well as their resolved URLs.
        template_example = url.replace("{release}", "110").replace("{identifier}", "ENST1")
        template_example = template_example.replace("{species}", "hsa").replace("{name}", "Mir-1")
        assert safe_source_url(template_example), url
    for forbidden in ("fetch(", "XMLHttpRequest", "localStorage", "sessionStorage", "document.cookie", "@import"):
        assert forbidden not in html, forbidden
