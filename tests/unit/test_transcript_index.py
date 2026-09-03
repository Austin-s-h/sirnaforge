"""Unit tests for the per-species transcript→gene index."""

from pathlib import Path

import pytest

from sirnaforge.data.transcript_index import (
    TranscriptGeneIndex,
    TranscriptRecord,
    _strip_version,
)


@pytest.fixture
def multi_species_fasta(tmp_path: Path) -> dict[str, Path]:
    """Create synthetic FASTA files for multiple species.

    Returns a dict mapping species names to their FASTA paths.
    """
    # Human FASTA with TP53 and a synthetic gene
    human_fasta = tmp_path / "human.cdna.fa"
    human_fasta.write_text(
        ">ENST00000269305.9 cdna chromosome:GRCh38:17:... "
        "gene:ENSG00000141510.16 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:TP53 "
        "description:tumor protein p53 [Source:HGNC Symbol;Acc:HGNC:11998]\n"
        "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTT\n"
        ">ENST00000987654.2 cdna chromosome:GRCh38:1:... "
        "gene:ENSG00000123456.7 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:SYNTH1 "
        "description:synthetic test gene 1\n"
        "ATGCGATCGATCGATCG\n"
        ">ENST00000111111.1 cdna chromosome:GRCh38:2:... "
        "gene:ENSG00000222222.3 gene_biotype:lncRNA "
        "transcript_biotype:lncRNA gene_symbol:LNCRNA1 "
        "description:long non-coding RNA 1\n"
        "GGGGCCCCAAAATTTT\n"
    )

    # Mouse FASTA with Trp53 (title case) and synthetic gene
    mouse_fasta = tmp_path / "mouse.cdna.fa"
    mouse_fasta.write_text(
        ">ENSMUST00000005493.13 cdna chromosome:GRCm39:11:... "
        "gene:ENSMUSG00000059552.17 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:Trp53 "
        "description:transformation related protein 53 [Source:MGI Symbol;Acc:MGI:98834]\n"
        "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTT\n"
        ">ENSMUST00000876543.4 cdna chromosome:GRCm39:5:... "
        "gene:ENSMUSG00000654321.2 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:SynthM1 "
        "description:synthetic mouse gene 1\n"
        "CCCCGGGGAAAATTTT\n"
        ">ENSMUST00000333333.1 cdna chromosome:GRCm39:7:... "
        "gene:ENSMUSG00000444444.5 gene_biotype:IG_V_gene "
        "gene_symbol:Igkv1-135 "
        "description:immunoglobulin kappa variable 1-135\n"
        "TTTTAAAACCCCGGGG\n"
    )

    # Rat FASTA with Tp53 and a transcript missing gene_symbol
    rat_fasta = tmp_path / "rat.cdna.fa"
    rat_fasta.write_text(
        ">ENSRNOT00000023349.7 cdna chromosome:mRatBN7.2:10:... "
        "gene:ENSRNOG00000017168.8 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:Tp53 "
        "description:tumor protein p53 [Source:RGD Symbol;Acc:RGD:3889]\n"
        "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTT\n"
        ">ENSRNOT00000999999.1 cdna chromosome:mRatBN7.2:3:... "
        "gene:ENSRNOG00000888888.2 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding "
        "description:gene without symbol annotation\n"
        "AAAACCCCGGGGTTTT\n"
    )

    return {
        "human": human_fasta,
        "mouse": mouse_fasta,
        "rat": rat_fasta,
    }


@pytest.mark.unit
def test_strip_version():
    """Version suffix stripping utility."""
    assert _strip_version("ENST00000269305.9") == "ENST00000269305"
    assert _strip_version("ENSG00000141510.16") == "ENSG00000141510"
    assert _strip_version("ENST00000269305") == "ENST00000269305"
    assert _strip_version("  ENST00000269305.9  ") == "ENST00000269305"


@pytest.mark.unit
def test_multi_species_parsing(multi_species_fasta: dict[str, Path]):
    """Parse gene IDs, symbols, and biotypes from multiple species."""
    index = TranscriptGeneIndex()

    # Build indices for all three species
    human_idx = index.build("human", multi_species_fasta["human"])
    mouse_idx = index.build("mouse", multi_species_fasta["mouse"])
    rat_idx = index.build("rat", multi_species_fasta["rat"])

    # Verify human parsing
    assert human_idx.species == "human"
    assert human_idx.transcript_count == 3
    assert human_idx.missing_symbol_count == 0

    tp53_record = human_idx.records.get("ENST00000269305")
    assert tp53_record is not None
    assert tp53_record.gene_id == "ENSG00000141510"
    assert tp53_record.gene_symbol == "TP53"
    assert tp53_record.biotype == "protein_coding"

    synth1_record = human_idx.records.get("ENST00000987654")
    assert synth1_record is not None
    assert synth1_record.gene_symbol == "SYNTH1"

    lnc_record = human_idx.records.get("ENST00000111111")
    assert lnc_record is not None
    assert lnc_record.biotype == "lncRNA"

    # Verify mouse parsing with title-case symbol
    assert mouse_idx.species == "mouse"
    assert mouse_idx.transcript_count == 3

    trp53_record = mouse_idx.records.get("ENSMUST00000005493")
    assert trp53_record is not None
    assert trp53_record.gene_id == "ENSMUSG00000059552"
    assert trp53_record.gene_symbol == "TRP53"  # Uppercased from Trp53

    igkv_record = mouse_idx.records.get("ENSMUST00000333333")
    assert igkv_record is not None
    assert igkv_record.biotype == "IG_V_gene"
    # transcript_biotype not present, falls back to gene_biotype via parse logic

    # Verify rat parsing
    assert rat_idx.species == "rat"
    assert rat_idx.transcript_count == 2
    assert rat_idx.missing_symbol_count == 1  # One record has no gene_symbol

    tp53_rat_record = rat_idx.records.get("ENSRNOT00000023349")
    assert tp53_rat_record is not None
    assert tp53_rat_record.gene_symbol == "TP53"  # Uppercased from Tp53


@pytest.mark.unit
def test_species_indices_do_not_collide(multi_species_fasta: dict[str, Path]):
    """Two species with overlapping transcript ID strings must not collide.

    This is the highest-value assertion: it prevents the regression where
    multiple species shared one dict.
    """
    # Add a collision scenario: same transcript ID prefix in both files
    tmp_path = multi_species_fasta["human"].parent
    collision_fasta = tmp_path / "collision_test.fa"
    collision_fasta.write_text(
        ">ENST00000999999.1 cdna gene:ENSG00000111111.1 "
        "gene_symbol:HUMANX gene_biotype:protein_coding "
        "transcript_biotype:protein_coding description:human gene\n"
        "AAAA\n"
    )

    collision_mouse_fasta = tmp_path / "collision_mouse_test.fa"
    collision_mouse_fasta.write_text(
        ">ENST00000999999.2 cdna gene:ENSMUSG00000222222.2 "
        "gene_symbol:MouseX gene_biotype:protein_coding "
        "transcript_biotype:protein_coding description:mouse gene\n"
        "CCCC\n"
    )

    index_collision = TranscriptGeneIndex()
    human_collision = index_collision.build("human", collision_fasta)
    mouse_collision = index_collision.build("mouse", collision_mouse_fasta)

    # The same transcript ID string "ENST00000999999" exists in both,
    # but each species' index resolves it to its own gene/symbol
    human_record = human_collision.records.get("ENST00000999999")
    mouse_record = mouse_collision.records.get("ENST00000999999")

    assert human_record is not None
    assert mouse_record is not None

    assert human_record.gene_id == "ENSG00000111111"
    assert human_record.gene_symbol == "HUMANX"

    assert mouse_record.gene_id == "ENSMUSG00000222222"
    assert mouse_record.gene_symbol == "MOUSEX"  # Uppercased

    # for_species returns the correct index
    retrieved_human = index_collision.for_species("human")
    assert retrieved_human is human_collision
    assert retrieved_human.symbol_for("ENST00000999999") == "HUMANX"

    retrieved_mouse = index_collision.for_species("mouse")
    assert retrieved_mouse is mouse_collision
    assert retrieved_mouse.symbol_for("ENST00000999999") == "MOUSEX"


@pytest.mark.unit
def test_version_suffix_tolerance(multi_species_fasta: dict[str, Path]):
    """Transcripts with version suffixes in FASTA are retrievable without suffix."""
    index = TranscriptGeneIndex()
    human_idx = index.build("human", multi_species_fasta["human"])

    # FASTA has ENST00000269305.9, but we query without version
    assert human_idx.gene_id_for("ENST00000269305") == "ENSG00000141510"
    assert human_idx.symbol_for("ENST00000269305") == "TP53"

    # Also works when querying with a version (it gets stripped)
    assert human_idx.gene_id_for("ENST00000269305.9") == "ENSG00000141510"
    assert human_idx.symbol_for("ENST00000269305.9") == "TP53"


@pytest.mark.unit
def test_case_insensitive_symbol_lookup(multi_species_fasta: dict[str, Path]):
    """Symbol lookup is case-insensitive across species."""
    index = TranscriptGeneIndex()
    human_idx = index.build("human", multi_species_fasta["human"])
    mouse_idx = index.build("mouse", multi_species_fasta["mouse"])
    rat_idx = index.build("rat", multi_species_fasta["rat"])

    # Human TP53 (all caps in FASTA)
    human_transcripts = human_idx.transcripts_for_symbol("tp53")
    assert "ENST00000269305" in human_transcripts

    # Also works with varying cases
    assert human_idx.transcripts_for_symbol("TP53") == human_transcripts
    assert human_idx.transcripts_for_symbol("Tp53") == human_transcripts

    # Mouse Trp53 (title case in FASTA, stored as TRP53)
    mouse_transcripts = mouse_idx.transcripts_for_symbol("trp53")
    assert "ENSMUST00000005493" in mouse_transcripts

    # Rat Tp53 (title case in FASTA, stored as TP53)
    rat_transcripts = rat_idx.transcripts_for_symbol("tp53")
    assert "ENSRNOT00000023349" in rat_transcripts

    # Cross-species comparison: all three species have a p53 ortholog
    # and they're all found via case-insensitive lookup
    assert len(human_transcripts) >= 1
    assert len(mouse_transcripts) >= 1
    assert len(rat_transcripts) >= 1


@pytest.mark.unit
def test_missing_gene_symbol_handling(multi_species_fasta: dict[str, Path]):
    """Records without gene_symbol return None and increment missing_symbol_count."""
    index = TranscriptGeneIndex()
    rat_idx = index.build("rat", multi_species_fasta["rat"])

    # Rat FASTA has one transcript without gene_symbol
    assert rat_idx.missing_symbol_count == 1

    # The transcript without a symbol
    no_symbol_record = rat_idx.records.get("ENSRNOT00000999999")
    assert no_symbol_record is not None
    assert no_symbol_record.gene_symbol is None
    assert rat_idx.symbol_for("ENSRNOT00000999999") is None

    # But it still has a gene_id and biotype
    assert no_symbol_record.gene_id == "ENSRNOG00000888888"
    assert no_symbol_record.biotype == "protein_coding"


@pytest.mark.unit
def test_description_with_spaces_and_colons(tmp_path: Path):
    """Headers where description: contains spaces and colons parse correctly."""
    fasta = tmp_path / "complex_desc.fa"
    fasta.write_text(
        ">ENST00000456789.3 cdna chromosome:GRCh38:7:... "
        "gene:ENSG00000112233.4 gene_biotype:protein_coding "
        "transcript_biotype:protein_coding gene_symbol:COMPLEX "
        "description:protein with complex name: domain A, variant B [Source:HGNC Symbol;Acc:HGNC:12345]\n"
        "ATGC\n"
    )

    index = TranscriptGeneIndex()
    idx = index.build("human", fasta)

    record = idx.records.get("ENST00000456789")
    assert record is not None
    assert record.gene_id == "ENSG00000112233"
    assert record.gene_symbol == "COMPLEX"
    assert record.biotype == "protein_coding"

    # The description field itself is not stored in TranscriptRecord,
    # but the earlier keys (gene, gene_symbol, etc.) must be parsed correctly
    # despite the description containing colons


@pytest.mark.unit
def test_absent_fasta_returns_empty_index(tmp_path: Path):
    """Missing FASTA path returns empty index without raising."""
    index = TranscriptGeneIndex()
    nonexistent = tmp_path / "does_not_exist.fa"

    idx = index.build("human", nonexistent)

    assert idx.species == "human"
    assert idx.transcript_count == 0
    assert idx.missing_symbol_count == 0
    assert len(idx.records) == 0

    # Lookups return None/empty
    assert idx.gene_id_for("ENST00000000000") is None
    assert idx.symbol_for("ENST00000000000") is None
    assert idx.transcripts_for_symbol("NOTFOUND") == frozenset()


@pytest.mark.unit
def test_idempotent_double_build(multi_species_fasta: dict[str, Path]):
    """Calling build twice for the same species returns the existing index."""
    index = TranscriptGeneIndex()

    # First build
    first = index.build("human", multi_species_fasta["human"])
    assert first.transcript_count == 3

    # Second build with the same species
    second = index.build("human", multi_species_fasta["human"])

    # Should return the same object
    assert second is first
    assert second.transcript_count == 3

    # No double-counting
    assert len(index.species()) == 1


@pytest.mark.unit
def test_index_lookup_methods(multi_species_fasta: dict[str, Path]):
    """Test all SpeciesTranscriptIndex lookup methods."""
    index = TranscriptGeneIndex()
    human_idx = index.build("human", multi_species_fasta["human"])

    # gene_id_for
    assert human_idx.gene_id_for("ENST00000269305") == "ENSG00000141510"
    assert human_idx.gene_id_for("NONEXISTENT") is None

    # symbol_for
    assert human_idx.symbol_for("ENST00000269305") == "TP53"
    assert human_idx.symbol_for("NONEXISTENT") is None

    # biotype_for
    assert human_idx.biotype_for("ENST00000269305") == "protein_coding"
    assert human_idx.biotype_for("ENST00000111111") == "lncRNA"
    assert human_idx.biotype_for("NONEXISTENT") is None

    # transcripts_for_symbol
    tp53_transcripts = human_idx.transcripts_for_symbol("TP53")
    assert "ENST00000269305" in tp53_transcripts
    assert len(tp53_transcripts) == 1

    synth1_transcripts = human_idx.transcripts_for_symbol("SYNTH1")
    assert "ENST00000987654" in synth1_transcripts

    # Non-existent symbol
    assert human_idx.transcripts_for_symbol("NOTREAL") == frozenset()


@pytest.mark.unit
def test_species_normalization(multi_species_fasta: dict[str, Path]):
    """Species names are normalized to canonical form."""
    index = TranscriptGeneIndex()

    # Build with various aliases
    idx1 = index.build("Homo sapiens", multi_species_fasta["human"])
    assert idx1.species == "human"

    # Retrieve with different alias
    idx2 = index.for_species("hsa")
    assert idx2 is idx1

    # Check species list
    assert "human" in index.species()


@pytest.mark.unit
def test_biotype_prefers_transcript_over_gene(tmp_path: Path):
    """Biotype prefers transcript_biotype, falling back to gene_biotype."""
    fasta = tmp_path / "biotype_test.fa"
    fasta.write_text(
        # Case 1: both present, prefer transcript_biotype
        ">ENST00000111111.1 cdna gene:ENSG00000111111.1 "
        "gene_biotype:protein_coding transcript_biotype:nonsense_mediated_decay "
        "gene_symbol:GENE1 description:test\n"
        "ATGC\n"
        # Case 2: only gene_biotype present
        ">ENST00000222222.1 cdna gene:ENSG00000222222.1 "
        "gene_biotype:lncRNA gene_symbol:GENE2 description:test\n"
        "ATGC\n"
    )

    index = TranscriptGeneIndex()
    idx = index.build("human", fasta)

    # Case 1: transcript_biotype wins
    record1 = idx.records.get("ENST00000111111")
    assert record1 is not None
    assert record1.biotype == "nonsense_mediated_decay"

    # Case 2: falls back to gene_biotype
    record2 = idx.records.get("ENST00000222222")
    assert record2 is not None
    assert record2.biotype == "lncRNA"


@pytest.mark.unit
def test_transcript_record_immutability():
    """TranscriptRecord is frozen and immutable."""
    record = TranscriptRecord(
        transcript_id="ENST00000269305",
        gene_id="ENSG00000141510",
        gene_symbol="TP53",
        biotype="protein_coding",
    )

    # Should not be able to modify (frozen dataclass)
    with pytest.raises((AttributeError, TypeError)):
        record.gene_symbol = "CHANGED"  # type: ignore
