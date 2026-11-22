"""Synthetic official-schema fixtures exercise public adapters without live services.

Schema references (these payloads are constructed, not recorded live responses):
https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/
https://rest.ensembl.org/documentation/info/variation_id
https://rest.ensembl.org/documentation/info/vep_hgvs_get
https://github.com/ncbi/dbsnp/blob/master/specs/refsnp_specification_deprecated.yaml
"""

import copy

import pytest

from sirnaforge.data.variant_resolution_helpers import VariantResolutionAmbiguity
from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import ClinVarSignificance, VariantSource

pytestmark = [pytest.mark.unit, pytest.mark.asyncio]


def _wire(resolver, monkeypatch, responses):
    """Replace only transport, leaving resolution and parsing code under test."""
    queue = iter(responses)
    calls = []

    async def request(url, params=None):
        calls.append((url, params))
        try:
            return copy.deepcopy(next(queue))
        except StopIteration:
            pytest.fail(f"Unexpected API request: {url}")

    monkeypatch.setattr(resolver, "_request_json", request)
    return calls


def _clinvar(uid, alt="A", assembly="GRCh38", accession="NC_000017.11", pos=100):
    """Construct ClinVar summary fields used to prove a genomic allele."""
    return {
        "uid": str(uid),
        "germline_classification": {"description": "Pathogenic"},
        "variation_set": [
            {
                "variation_xrefs": [{"db_source": "dbSNP", "db_id": "12"}],
                "canonical_spdi": f"{accession}:{pos - 1}:G:{alt}",
                "variation_loc": [
                    {"assembly_name": assembly, "chr": "17", "start": str(pos), "stop": str(pos), "ref": "", "alt": ""}
                ],
            }
        ],
    }


def _ensembl(alleles="G/A", **mapping):
    """Construct documented Ensembl mapping and allele population fields."""
    return {
        "name": "rs12",
        "mappings": [
            {
                "assembly_name": "GRCh38",
                "seq_region_name": "17",
                "start": 100,
                "end": 100,
                "strand": 1,
                "allele_string": alleles,
                **mapping,
            }
        ],
        "populations": [
            {"population": "gnomADe:ALL", "allele": "G", "frequency": 0.8},
            {"population": "gnomADe:ALL", "allele": "A", "frequency": 0.12},
            {"population": "gnomADe:ALL", "allele": "T", "frequency": 0.08},
            {"population": "gnomADe:AFR", "allele": "T", "frequency": 0.3},
        ],
    }


def _vep(hgvs, pos=100, **extra):
    """Construct HGVS VEP output including explicit GRCh38 and VCF alleles."""
    return {
        "input": hgvs,
        "assembly_name": "GRCh38",
        "seq_region_name": "17",
        "start": pos,
        "end": pos,
        "strand": 1,
        "allele_string": "G/A",
        "vcf_string": f"17-{pos}-G-A",
        **extra,
    }


def _refsnp(alt="A"):
    """Construct RefSNP SPDI placements with parallel allele annotations."""
    placement = {
        "seq_id": "NC_000017.11",
        "placement_annot": {
            "seq_id_traits_by_assembly": [
                {
                    "assembly_name": "GRCh38.p14",
                    "is_top_level": True,
                    "is_chromosome": True,
                    "is_alt": False,
                    "is_patch": False,
                }
            ]
        },
        "alleles": [
            {
                "allele": {
                    "spdi": {
                        "seq_id": "NC_000017.11",
                        "position": 99,
                        "deleted_sequence": "G",
                        "inserted_sequence": allele,
                    }
                }
            }
            for allele in ("G", alt)
        ],
    }
    return {
        "refsnp_id": "12",
        "primary_snapshot_data": {
            "placements_with_allele": [placement],
            "allele_annotations": [
                {"frequency": [{"study_name": "GnomAD", "allele_count": count, "total_count": 100}]}
                for count in (90, 10)
            ],
        },
    }


async def test_clinvar_uses_chrpos38_and_exact_second_hit(tmp_path, monkeypatch):
    """A nearby first hit cannot replace the requested exact alternate allele."""
    resolver = VariantResolver(cache_dir=tmp_path)
    calls = _wire(
        resolver,
        monkeypatch,
        [
            {"esearchresult": {"idlist": ["1", "2"], "count": "2"}},
            {"result": {"1": _clinvar(1, alt="T"), "2": _clinvar(2)}},
        ],
    )
    result = await resolver._query_clinvar(resolver.parse_identifier("17:100:G:A"))
    assert result is not None
    assert result.to_vcf_style() == "chr17:100:G:A"
    assert result.provenance["clinvar_ids"] == ["2"]
    assert calls[0][1]["term"] == "17[chr] AND 100[chrpos38]"
    assert calls[1][1]["id"] == "1,2"


@pytest.mark.parametrize(
    "record",
    [_clinvar(1, assembly="GRCh37"), _clinvar(1, accession="NC_000017.10"), _clinvar(1, pos=101), _clinvar(1, alt="T")],
)
async def test_clinvar_rejects_unproven_assembly_locus_and_alleles(tmp_path, monkeypatch, record):
    """Returned coordinates and alleles must independently match the query."""
    resolver = VariantResolver(cache_dir=tmp_path)
    _wire(resolver, monkeypatch, [{"esearchresult": {"idlist": ["1"]}}, {"result": {"1": record}}])
    assert await resolver._query_clinvar(resolver.parse_identifier("17:100:G:A")) is None


async def test_clinvar_does_not_invent_missing_alleles_from_input(tmp_path, monkeypatch):
    """A matching location without actual alleles is insufficient evidence."""
    resolver = VariantResolver(cache_dir=tmp_path)
    record = _clinvar(1)
    del record["variation_set"][0]["canonical_spdi"]
    _wire(resolver, monkeypatch, [{"esearchresult": {"idlist": ["1"]}}, {"result": {"1": record}}])
    assert await resolver._query_clinvar(resolver.parse_identifier("17:100:G:A")) is None


async def test_clinvar_ambiguous_rsid_stops_source_fallback(tmp_path, monkeypatch):
    """An ambiguous rsID requires an explicit allele instead of a lower-source guess."""
    resolver = VariantResolver(cache_dir=tmp_path)
    calls = _wire(
        resolver,
        monkeypatch,
        [{"esearchresult": {"idlist": ["1", "2"]}}, {"result": {"1": _clinvar(1), "2": _clinvar(2, alt="T")}}],
    )
    assert await resolver.resolve_variant(resolver.parse_identifier("rs12")) is None
    assert len(calls) == 2


async def test_clinvar_truncated_search_cannot_choose_partial_hit_list(tmp_path, monkeypatch):
    """A bounded query must identify truncation instead of silently choosing a subset."""
    resolver = VariantResolver(cache_dir=tmp_path)
    _wire(resolver, monkeypatch, [{"esearchresult": {"idlist": ["1"], "count": "101"}}])
    with pytest.raises(VariantResolutionAmbiguity, match="100-hit"):
        await resolver._query_clinvar(resolver.parse_identifier("rs12"))


async def test_clinvar_conflicting_duplicate_classifications_fail_default_filter(tmp_path, monkeypatch):
    """Matching alleles retain duplicate provenance and expose classification conflicts."""
    resolver = VariantResolver(cache_dir=tmp_path)
    second = _clinvar(2)
    second["germline_classification"]["description"] = "Benign"
    _wire(
        resolver, monkeypatch, [{"esearchresult": {"idlist": ["1", "2"]}}, {"result": {"1": _clinvar(1), "2": second}}]
    )
    result = await resolver._query_clinvar(resolver.parse_identifier("17:100:G:A"))
    assert result is not None
    assert result.clinvar_significance == ClinVarSignificance.CONFLICTING
    assert result.provenance["clinvar_ids"] == ["1", "2"]
    assert resolver._passes_filters(result) is False


async def test_ensembl_skips_wrong_assembly_mapping_and_reads_alt_specific_af(tmp_path, monkeypatch):
    """Neither a first GRCh37 mapping nor reference-allele frequency can leak through."""
    resolver = VariantResolver(cache_dir=tmp_path)
    data = _ensembl()
    data["mappings"].insert(0, {**data["mappings"][0], "assembly_name": "GRCh37", "start": 90, "end": 90})
    _wire(resolver, monkeypatch, [data])
    result = await resolver._query_ensembl(resolver.parse_identifier("rs12"))
    assert result is not None
    assert result.to_vcf_style() == "chr17:100:G:A"
    assert result.af == 0.12


async def test_ensembl_multiallelic_rsid_requires_explicit_alternate(tmp_path, monkeypatch):
    """Single-record API cannot choose the first ALT from a multiallelic rsID."""
    resolver = VariantResolver(cache_dir=tmp_path)
    _wire(resolver, monkeypatch, [_ensembl("G/A/T")])
    with pytest.raises(VariantResolutionAmbiguity):
        await resolver._query_ensembl(resolver.parse_identifier("rs12"))


async def test_ensembl_overlap_uses_allele_list_and_enriches_selected_alt_only(tmp_path, monkeypatch):
    """Real overlap allele arrays support exact selection and allele-specific frequencies."""
    resolver = VariantResolver(cache_dir=tmp_path)
    overlap = [
        {
            "id": "rs12",
            "assembly_name": "GRCh38",
            "seq_region_name": "17",
            "start": 100,
            "end": 100,
            "strand": 1,
            "alleles": ["G", "A", "T"],
        }
    ]
    _wire(resolver, monkeypatch, [overlap, _ensembl("G/A/T")])
    result = await resolver._query_ensembl(resolver.parse_identifier("17:100:G:T"))
    assert result is not None
    assert result.alt == "T"
    assert result.af == 0.08
    assert result.population_afs["AFR"] == 0.3


async def test_ensembl_negative_strand_returns_forward_alleles(tmp_path, monkeypatch):
    """Alleles are reverse-complemented together while frequencies stay attached to source ALT."""
    resolver = VariantResolver(cache_dir=tmp_path)
    _wire(resolver, monkeypatch, [_ensembl(strand=-1)])
    result = await resolver._query_ensembl(resolver.parse_identifier("rs12"))
    assert result is not None
    assert result.to_vcf_style() == "chr17:100:C:T"
    assert result.af == 0.12


async def test_hgvs_normalizes_before_clinvar_priority_and_reuses_cache(tmp_path, monkeypatch):
    """HGVS resolution produces an exact allele and preserves ClinVar as the selected source."""
    resolver = VariantResolver(cache_dir=tmp_path)
    hgvs = "NM_000546.6:c.215C>G"
    calls = _wire(
        resolver, monkeypatch, [[_vep(hgvs)], {"esearchresult": {"idlist": ["1"]}}, {"result": {"1": _clinvar(1)}}]
    )
    query = resolver.parse_identifier(hgvs)
    result = await resolver.resolve_variant(query)
    assert result is not None
    assert result.to_vcf_style() == "chr17:100:G:A"
    assert result.sources == [VariantSource.CLINVAR]
    assert result.provenance["hgvs_normalization"]["input"] == hgvs
    assert calls[0][1] == {"vcf_string": "1", "ambiguous_hgvs": "1"}
    assert calls[1][1]["term"] == "17[chr] AND 100[chrpos38]"
    assert await resolver.resolve_variant(query) == result
    assert len(calls) == 3


async def test_hgvs_unregistered_variant_can_use_verified_vep_allele(tmp_path, monkeypatch):
    """A normalized allele need not already have an rsID to support HGVS input."""
    resolver = VariantResolver(cache_dir=tmp_path, source_priority=[VariantSource.ENSEMBL])
    hgvs = "NM_000546.6:c.215C>G"
    _wire(resolver, monkeypatch, [[_vep(hgvs)], []])
    result = await resolver.resolve_variant(resolver.parse_identifier(hgvs))
    assert result is not None
    assert result.id == hgvs
    assert result.sources == [VariantSource.ENSEMBL]


async def test_ambiguous_protein_hgvs_does_not_pick_first_mapping(tmp_path, monkeypatch):
    """Protein descriptions mapping to different nucleotide alleles remain unresolved."""
    resolver = VariantResolver(cache_dir=tmp_path)
    hgvs = "NP_000537.3:p.Pro72Arg"
    calls = _wire(resolver, monkeypatch, [[_vep(hgvs), _vep(hgvs, pos=200)]])
    assert await resolver.resolve_variant(resolver.parse_identifier(hgvs)) is None
    assert len(calls) == 1


@pytest.mark.parametrize(
    "override",
    [
        {"assembly_name": "GRCh37"},
        {"input": "different"},
        {"seq_region_name": "GL000009.2"},
        {"vcf_string": "17-100-G-A,T"},
    ],
)
async def test_hgvs_rejects_unproven_or_nonprimary_normalization(tmp_path, monkeypatch, override):
    """HGVS mappings require exact input identity, assembly and a single primary allele."""
    resolver = VariantResolver(cache_dir=tmp_path)
    hgvs = "NM_000546.6:c.215C>G"
    _wire(resolver, monkeypatch, [[_vep(hgvs, **override)]])
    assert await resolver.resolve_variant(resolver.parse_identifier(hgvs)) is None


async def test_dbsnp_rsid_reads_grch38_spdi_and_alt_annotation_frequency(tmp_path, monkeypatch):
    """RefSNP offsets and per-allele parallel annotations produce the correct record."""
    resolver = VariantResolver(cache_dir=tmp_path)
    data = _refsnp()
    wrong = copy.deepcopy(data["primary_snapshot_data"]["placements_with_allele"][0])
    wrong["seq_id"] = "NC_000017.10"
    wrong["placement_annot"]["seq_id_traits_by_assembly"][0]["assembly_name"] = "GRCh37"
    data["primary_snapshot_data"]["placements_with_allele"].insert(0, wrong)
    _wire(resolver, monkeypatch, [data])
    result = await resolver._query_dbsnp(resolver.parse_identifier("rs12"))
    assert result is not None
    assert result.to_vcf_style() == "chr17:100:G:A"
    assert result.af == 0.1
    assert result.provenance["frequency_source"] == "GnomAD"


async def test_dbsnp_coordinate_maps_exact_spdi_before_fetch(tmp_path, monkeypatch):
    """Coordinate fallback uses the pinned GRCh38 chromosome accession and exact allele."""
    resolver = VariantResolver(cache_dir=tmp_path)
    calls = _wire(resolver, monkeypatch, [{"data": {"rsids": [12]}}, _refsnp()])
    result = await resolver._query_dbsnp(resolver.parse_identifier("17:100:G:A"))
    assert result is not None
    assert "NC_000017.11%3A99%3AG%3AA/rsids" in calls[0][0]


async def test_dbsnp_does_not_accept_neighbour_from_exact_lookup(tmp_path, monkeypatch):
    """The returned RefSNP placement still must match even after exact-SPDI lookup."""
    resolver = VariantResolver(cache_dir=tmp_path)
    _wire(resolver, monkeypatch, [{"data": {"rsids": [12]}}, _refsnp(alt="T")])
    assert await resolver._query_dbsnp(resolver.parse_identifier("17:100:G:A")) is None


async def test_default_source_chain_reaches_real_dbsnp_adapter(tmp_path, monkeypatch):
    """Missing ClinVar and Ensembl evidence can fall back to a verified RefSNP allele."""
    resolver = VariantResolver(cache_dir=tmp_path)
    calls = _wire(resolver, monkeypatch, [{"esearchresult": {"idlist": []}}, {}, _refsnp()])
    result = await resolver.resolve_variant(resolver.parse_identifier("rs12"))
    assert result is not None
    assert result.sources == [VariantSource.DBSNP]
    assert len(calls) == 3


async def test_source_priority_order_partitions_cache(tmp_path):
    """Changing trusted source order must never serve the previous order's cached winner."""
    first = VariantResolver(cache_dir=tmp_path, source_priority=[VariantSource.CLINVAR, VariantSource.ENSEMBL])
    second = VariantResolver(cache_dir=tmp_path, source_priority=[VariantSource.ENSEMBL, VariantSource.CLINVAR])
    query = first.parse_identifier("rs12")
    assert first._get_cache_key(query) != second._get_cache_key(query)


async def test_combined_pathogenic_classification_is_accepted_without_losing_raw_label(tmp_path, monkeypatch):
    """ClinVar's combined pathogenic label meets the default two-category policy."""
    resolver = VariantResolver(cache_dir=tmp_path)
    item = _clinvar(1)
    item["germline_classification"]["description"] = "Pathogenic/Likely pathogenic"
    _wire(resolver, monkeypatch, [{"esearchresult": {"idlist": ["1"]}}, {"result": {"1": item}}])
    result = await resolver.resolve_variant(resolver.parse_identifier("17:100:G:A"))
    assert result is not None
    assert result.clinvar_significance.value == "Pathogenic/Likely pathogenic"
