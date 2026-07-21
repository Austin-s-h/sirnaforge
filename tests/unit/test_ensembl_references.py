"""Tests for the shared Ensembl reference table driving transcriptome + genome sources."""

from sirnaforge.data.ensembl_references import (
    ENSEMBL_ASSEMBLIES,
    build_genome_sources,
    build_transcriptome_sources,
)
from sirnaforge.data.genome_manager import GenomeManager
from sirnaforge.data.transcriptome_manager import TranscriptomeManager

EXPECTED_SPECIES = {"human", "mouse", "rat", "macaque"}


def test_table_covers_expected_species() -> None:
    """The assembly table should cover the same species set as the transcriptome/miRNA layers."""
    assert {assembly.species for assembly in ENSEMBL_ASSEMBLIES} == EXPECTED_SPECIES


def test_transcriptome_sources_preserve_legacy_keys_and_urls() -> None:
    """Generated cDNA sources must keep the historical keys, names, URLs, and cache keys.

    The cache key hashes name+species+url, so any drift would silently invalidate every
    user's existing transcriptome cache. These URLs are frozen intentionally.
    """
    sources = build_transcriptome_sources()
    assert set(sources) == {f"ensembl_{sp}_cdna" for sp in EXPECTED_SPECIES}

    expected_urls = {
        "ensembl_human_cdna": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "ensembl_mouse_cdna": "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
        "ensembl_rat_cdna": "https://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz",
        "ensembl_macaque_cdna": "https://ftp.ensembl.org/pub/current_fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz",
    }
    for key, url in expected_urls.items():
        assert sources[key].url == url
        assert sources[key].name == "ensembl_cdna"


def test_transcriptome_manager_uses_generated_sources() -> None:
    """TranscriptomeManager.SOURCES should equal the generated table (same keys)."""
    assert set(TranscriptomeManager.SOURCES) == set(build_transcriptome_sources())


def test_genome_sources_expose_all_species_and_preserve_human_key() -> None:
    """Genome sources should cover all species and keep the legacy human key."""
    sources = build_genome_sources()
    assert set(sources) == {
        "ensembl_human_hg38_primary",
        "ensembl_mouse_grcm39_primary",
        "ensembl_rat_grcr8_toplevel",
        "ensembl_macaque_mmul10_toplevel",
    }
    # Backward-compatibility: the ZFN default / CLI / existing caches reference this key.
    assert "ensembl_human_hg38_primary" in sources
    assert sources["ensembl_human_hg38_primary"].species == "human"


def test_genome_dna_urls_use_correct_per_species_file_kind() -> None:
    """Human/mouse publish primary_assembly; rat/macaque only publish toplevel.

    A wrong dna_kind is a silent 404 at download time, so this pins the verified layout.
    """
    sources = build_genome_sources()
    assert "dna.primary_assembly.fa.gz" in sources["ensembl_human_hg38_primary"].url
    assert "dna.primary_assembly.fa.gz" in sources["ensembl_mouse_grcm39_primary"].url
    assert "dna.toplevel.fa.gz" in sources["ensembl_rat_grcr8_toplevel"].url
    assert "dna.toplevel.fa.gz" in sources["ensembl_macaque_mmul10_toplevel"].url


def test_genome_manager_uses_generated_sources() -> None:
    """GenomeManager should surface the generated multi-species genome sources."""
    manager = GenomeManager(auto_build_indices=False)
    assert set(manager.sources) == set(build_genome_sources())


def test_cdna_and_genome_share_species_slug_and_assembly() -> None:
    """Both cDNA and DNA URLs for a species must share the same Ensembl slug + assembly."""
    for assembly in ENSEMBL_ASSEMBLIES:
        cdna = assembly.cdna_url()
        dna = assembly.dna_url()
        assert f"/{assembly.ensembl_slug}/" in cdna
        assert f"/{assembly.ensembl_slug}/" in dna
        assert assembly.assembly in cdna
        assert assembly.assembly in dna
