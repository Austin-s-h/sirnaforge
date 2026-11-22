"""Deterministic source links require stable identifiers and sanitized resolver templates."""

from __future__ import annotations

import pytest

from sirnaforge.reporting.links import canonical_identity_link, safe_source_url, source_link


@pytest.mark.unit
def test_ensembl_resolver_preserves_identifier_version_and_release() -> None:
    """Release aliases provide a version-aware link without a render-time network request."""
    link = source_link("ENST00000359597.8", "ensembl", release="110")
    assert link["url"] == "https://e110.ensembl.org/id/ENST00000359597.8"
    assert link["release"] == "110"
    current = source_link("ENSG00000141510", "ensembl")
    assert current["url"] == "https://www.ensembl.org/id/ENSG00000141510"
    assert current["release"] is None


@pytest.mark.unit
def test_manifest_release_template_takes_precedence_over_current_template() -> None:
    """A manifest may declare both aliases without losing a version-pinned Ensembl link."""
    resolver = {
        "template": "https://www.ensembl.org/id/{identifier}",
        "release_template": "https://e{release}.ensembl.org/id/{identifier}",
    }
    pinned = source_link("ENST1", "ensembl", release="110", resolver=resolver)
    current = source_link("ENST1", "ensembl", resolver=resolver)
    assert pinned["url"] == "https://e110.ensembl.org/id/ENST1"
    assert current["url"] == "https://www.ensembl.org/id/ENST1"


@pytest.mark.unit
def test_gencode_annotation_release_is_not_used_as_ensembl_release() -> None:
    """GENCODE's annotation version identifies the authority, not Ensembl's URL release."""
    transcript = canonical_identity_link("ENST1.2", "GENCODE explicit canonical/MANE tags (annotation 49)", "49")
    gene = canonical_identity_link("ENSG1.3", "GENCODE annotation 49", "49")
    assert transcript["release"] is None
    assert gene["release"] is None
    assert transcript["url"] == "https://www.ensembl.org/id/ENST1.2"
    assert gene["url"] == "https://www.ensembl.org/id/ENSG1.3"


@pytest.mark.unit
def test_named_mirgenedb_mature_hit_resolves_to_its_locus() -> None:
    """Named mature arm identifiers map to their actual locus, preserving the original ID."""
    link = source_link("Hsa-Mir-24-P2_3p", "mirgenedb")
    assert link["url"] == "https://mirgenedb.org/show/hsa/Mir-24-P2"
    assert link["identifier"] == "Hsa-Mir-24-P2_3p"
    assert source_link("hsa-miR-24-3p", "mirgenedb")["url"] is None
    assert source_link("TP53", "ensembl")["url"] is None


@pytest.mark.unit
@pytest.mark.parametrize(
    "template",
    [
        "javascript:alert(1)",
        "http://www.ensembl.org/id/{identifier}",
        "https://www.ensembl.org.evil.example/id/{identifier}",
        "https://evil.example/id/{identifier}",
        "https://user:pass@www.ensembl.org/id/{identifier}",
        'https://www.ensembl.org/id/{identifier}" onclick="alert(1)',
        "https://www.ensembl.org/id/{missing}",
    ],
)
def test_invalid_source_resolvers_do_not_produce_links(template: str) -> None:
    """A manifest resolver cannot introduce active content or an unapproved external destination."""
    assert source_link("ENST1", "ensembl", resolver={"template": template})["url"] is None
    assert not safe_source_url("https://www.ensembl.org:invalid/id/ENST1")
