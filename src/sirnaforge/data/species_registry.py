"""Canonical species registry and metadata for miRNA and genome mappings."""

from __future__ import annotations

from typing import Any

# MirGeneDB metadata built around canonical species slugs (e.g. hsa, mmu).
MIRGENEDB_SPECIES_TABLE: dict[str, dict[str, Any]] = {
    "hsa": {
        "taxonomy_id": "9606",
        "scientific_name": "Homo sapiens",
        "common_name": "Human",
        "aliases": ["human", "homo_sapiens"],
    },
    "ptr": {
        "taxonomy_id": "9598",
        "scientific_name": "Pan troglodytes",
        "common_name": "Chimpanzee",
        "aliases": ["chimpanzee", "chimp", "pan_troglodytes"],
    },
    "ggo": {
        "taxonomy_id": "9593",
        "scientific_name": "Gorilla gorilla",
        "common_name": "Gorilla",
        "aliases": ["gorilla", "gorilla_gorilla"],
    },
    "mml": {
        "taxonomy_id": "9544",
        "scientific_name": "Macaca mulatta",
        "common_name": "Rhesus macaque",
        "aliases": ["rhesus", "macaca_mulatta"],
    },
    "mmu": {
        "taxonomy_id": "10090",
        "scientific_name": "Mus musculus",
        "common_name": "Mouse",
        "aliases": ["mouse", "mus_musculus"],
    },
    "rno": {
        "taxonomy_id": "10116",
        "scientific_name": "Rattus norvegicus",
        "common_name": "Rat",
        "aliases": ["rat", "rattus_norvegicus"],
    },
    "cfa": {
        "taxonomy_id": "9615",
        "scientific_name": "Canis lupus familiaris",
        "common_name": "Dog",
        "aliases": ["dog", "canis_lupus_familiaris"],
    },
    "fca": {
        "taxonomy_id": "9685",
        "scientific_name": "Felis catus",
        "common_name": "Cat",
        "aliases": ["cat", "felis_catus"],
    },
    "bta": {
        "taxonomy_id": "9913",
        "scientific_name": "Bos taurus",
        "common_name": "Cow",
        "aliases": ["cow", "bos_taurus"],
    },
    "ssc": {
        "taxonomy_id": "9823",
        "scientific_name": "Sus scrofa",
        "common_name": "Pig",
        "aliases": ["pig", "sus_scrofa"],
    },
    "eca": {
        "taxonomy_id": "9796",
        "scientific_name": "Equus caballus",
        "common_name": "Horse",
        "aliases": ["horse", "equus_caballus"],
    },
    "oar": {
        "taxonomy_id": "9940",
        "scientific_name": "Ovis aries",
        "common_name": "Sheep",
        "aliases": ["sheep", "ovis_aries"],
    },
    "gga": {
        "taxonomy_id": "9031",
        "scientific_name": "Gallus gallus",
        "common_name": "Chicken",
        "aliases": ["chicken", "gallus_gallus", "chi"],
    },
    "tgu": {
        "taxonomy_id": "9103",
        "scientific_name": "Meleagris gallopavo",
        "common_name": "Turkey",
        "aliases": ["turkey", "meleagris_gallopavo"],
    },
    "dre": {
        "taxonomy_id": "7955",
        "scientific_name": "Danio rerio",
        "common_name": "Zebrafish",
        "aliases": ["zebrafish", "danio_rerio"],
    },
    "xla": {
        "taxonomy_id": "8355",
        "scientific_name": "Xenopus laevis",
        "common_name": "African clawed frog",
        "aliases": ["xenopus", "xenopus_laevis"],
    },
    "ola": {
        "taxonomy_id": "8090",
        "scientific_name": "Oryzias latipes",
        "common_name": "Medaka",
        "aliases": ["medaka", "oryzias_latipes"],
    },
    "gac": {
        "taxonomy_id": "69293",
        "scientific_name": "Gasterosteus aculeatus",
        "common_name": "Stickleback",
        "aliases": ["stickleback", "gasterosteus_aculeatus"],
    },
    "pma": {
        "taxonomy_id": "7757",
        "scientific_name": "Petromyzon marinus",
        "common_name": "Sea lamprey",
        "aliases": ["lamprey", "petromyzon_marinus"],
    },
    "dme": {
        "taxonomy_id": "7227",
        "scientific_name": "Drosophila melanogaster",
        "common_name": "Fruit fly",
        "aliases": ["fruitfly", "drosophila_melanogaster", "dmel"],
    },
    "aga": {
        "taxonomy_id": "7165",
        "scientific_name": "Anopheles gambiae",
        "common_name": "Mosquito",
        "aliases": ["mosquito", "anopheles_gambiae"],
    },
    "cel": {
        "taxonomy_id": "6239",
        "scientific_name": "Caenorhabditis elegans",
        "common_name": "C. elegans",
        "aliases": ["worm", "caenorhabditis_elegans"],
    },
    "spur": {
        "taxonomy_id": "7668",
        "scientific_name": "Strongylocentrotus purpuratus",
        "common_name": "Purple sea urchin",
        "aliases": ["seaurchin", "sea_urchin", "strongylocentrotus_purpuratus"],
    },
}

MIRGENEDB_ALIAS_MAP: dict[str, str] = {}
for slug, metadata in MIRGENEDB_SPECIES_TABLE.items():
    MIRGENEDB_ALIAS_MAP[slug.lower()] = slug
    for alias in metadata.get("aliases", []):
        MIRGENEDB_ALIAS_MAP[alias.lower()] = slug


CANONICAL_SPECIES_REGISTRY: dict[str, dict[str, Any]] = {
    "human": {
        "genome": "human",
        "mirgenedb_slug": "hsa",
        "aliases": ["homo_sapiens"],
    },
    "mouse": {
        "genome": "mouse",
        "mirgenedb_slug": "mmu",
        "aliases": ["mus_musculus"],
    },
    "rhesus": {
        "genome": "rhesus",
        "mirgenedb_slug": "mml",
        "aliases": ["macaca_mulatta", "rhesus_macaque"],
    },
    "rat": {
        "genome": "rat",
        "mirgenedb_slug": "rno",
        "aliases": ["rattus_norvegicus"],
    },
    "chicken": {
        "genome": "chicken",
        "mirgenedb_slug": "gga",
        "aliases": ["gallus_gallus", "chi"],
    },
}

CANONICAL_SPECIES_ALIAS_MAP: dict[str, str] = {}
for canonical_name, registry_entry in CANONICAL_SPECIES_REGISTRY.items():
    alias_values = {canonical_name}
    alias_values.update({alias.lower() for alias in registry_entry.get("aliases", [])})

    slug = registry_entry["mirgenedb_slug"]
    alias_values.add(slug.lower())

    mirgenedb_metadata = MIRGENEDB_SPECIES_TABLE.get(slug)
    if mirgenedb_metadata:
        alias_values.update({alias.lower() for alias in mirgenedb_metadata.get("aliases", [])})
        common_name = mirgenedb_metadata.get("common_name")
        if common_name:
            alias_values.add(common_name.lower())
        scientific_name = mirgenedb_metadata.get("scientific_name")
        if scientific_name:
            alias_values.add(scientific_name.lower())

    for alias in alias_values:
        if alias:
            CANONICAL_SPECIES_ALIAS_MAP[alias] = canonical_name
