# Modular Override Guide

## Design Philosophy

siRNAforge uses a **single canonical species selector** (`--species`) that automatically maps to all off-target analysis layers:

```bash
# All formats work - automatically normalized to canonical names
sirnaforge workflow TP53 --species human,mouse,rat
sirnaforge workflow TP53 --species hsa,mmu,rno           # miRBase codes
sirnaforge workflow TP53 --species "Homo sapiens,Mus musculus,Rattus norvegicus"
```

**Automatic Normalization:** Species names are intelligently converted to canonical form (`human`, `mouse`, `rat`, etc.) regardless of input format:
- **Common names**: `human`, `mouse`, `macaque`, `rhesus`
- **miRBase codes**: `hsa`, `mmu`, `rno`, `mml`, `gga`, `ssc`
- **Scientific names**: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`

This single parameter drives:
- **miRNA database lookups** → MirGeneDB/miRBase queries for all specified species
- **Transcriptome fetching** → Automatic Ensembl cDNA downloads for available species
- **Report aggregation** → Consolidated off-target summaries across species

You can **surgically override** any layer without affecting others.

---

## Default Behavior

### Without Any Overrides

```bash
sirnaforge workflow TP53
```

**What happens:**
1. **Species**: Uses `DEFAULT_MIRNA_CANONICAL_SPECIES` (7 species: chicken, pig, rat, mouse, human, rhesus, macaque)
2. **miRNA Database**: MirGeneDB lookups for all 7 species
3. **Transcriptome**: Ensembl cDNA auto-fetched for 4 species (human, mouse, rat, macaque)
4. **Genomic DNA**: None (optional, resource-intensive)

---

## Surgical Overrides

### Override miRNA Species Only

```bash
# Use human species everywhere, but only check rat miRNAs
sirnaforge workflow TP53 \
  --species human \
  --mirna-species rno
```

**Result:** Human transcriptome + rat miRNA database

### Override Transcriptome Reference

```bash
# Custom transcriptome for human analysis
sirnaforge workflow TP53 \
  --species human \
  --transcriptome-fasta /path/to/custom_transcripts.fa
```

**Options for `--transcriptome-fasta`:**
- **Pre-configured source**: `ensembl_human_cdna`, `ensembl_mouse_cdna`, etc.
- **Local file**: `/path/to/transcripts.fa`
- **Remote URL**: `https://ftp.ensembl.org/.../cdna.fa.gz`

The system **automatically**:
- Downloads remote files to cache
- Decompresses `.gz` archives
- Builds BWA-MEM2 indices
- Reuses indices on subsequent runs

### Add Custom Sequences to Defaults

**Current behavior**: `--transcriptome-fasta` *replaces* defaults.

**Additive capability** (implementation detail): The `transcriptome_selection.choices` tuple can hold multiple `ReferenceChoice` objects, allowing future support for comma-separated lists:

```bash
# Future enhancement:
sirnaforge workflow TP53 \
  --transcriptome-fasta ensembl_human_cdna,/path/to/custom_contig.fa
```

**Current workaround**: Manually concatenate FASTA files:

```bash
cat ensembl_human_cdna.fa custom_sequences.fa > combined.fa
sirnaforge workflow TP53 --transcriptome-fasta combined.fa
```

### Override BWA Indices Directly

Skip automatic fetching and use pre-built indices:

```bash
sirnaforge workflow TP53 \
  --offtarget-indices human:/data/GRCh38_index,mouse:/data/GRCm39_index
```

**Format**: `species:/absolute/path/to/index_prefix`
**Effect**: Bypasses Ensembl downloads and uses your existing BWA-MEM2 indices

---

## Advanced Use Cases

### Custom miRNA Database (Not Recommended)

The miRNA system doesn't currently support custom FASTA override via CLI, but you can:

1. **Use the API** directly with `MiRNADatabaseManager.get_custom_database()` (undocumented)
2. **Add sequences to cache** by placing files in `~/.cache/sirnaforge/mirna/`
3. **Request the feature** if this is a common workflow

**Why discouraged?** miRNA databases use species-specific three-letter codes (e.g., `hsa`, `mmu`) and specialized filtering. Using canonical sources (MirGeneDB, miRBase) ensures consistency.

### Design-Only Mode

Skip all off-target analysis:

```bash
sirnaforge design transcripts.fasta --top-n 50
```

**Result:** Generates `_all.csv` and `_pass.csv` with thermodynamic scores only.

### Genomic DNA Off-Target (Resource-Intensive)

Requires pre-built BWA indices for full genomes:

```bash
sirnaforge workflow TP53 \
  --offtarget-indices human:/data/genomes/GRCh38_bwa
```

**Storage requirements:** ~3GB per genome index
**Compute requirements:** ~4-8GB RAM per parallel BWA job

---

## Override Hierarchy

When multiple parameters affect the same resource:

1. **Explicit overrides win**: `--offtarget-indices` > `--transcriptome-fasta` > `--species`
2. **miRNA overrides are independent**: `--mirna-species` doesn't affect transcriptomes
3. **Defaults are smart**: System auto-detects what's available (e.g., only 4/7 species have Ensembl cDNA)

### Example: Mixed Override

```bash
sirnaforge workflow TP53 \
  --species human,mouse,rat,macaque \
  --mirna-species hsa,mmu,rno,mml,gga,ssc \
  --transcriptome-fasta /custom/human_isoforms.fa \
  --offtarget-indices mouse:/data/GRCm39,rat:/data/rn7
```

**Interpretation:**
- **Canonical species**: 4 (human, mouse, rat, macaque)
- **miRNA checks**: 6 species (includes chicken, pig)
- **Transcriptome**: Custom human file (replaces default Ensembl human cDNA)
- **Off-target indices**: Custom mouse/rat indices (replaces Ensembl-fetched)
- **Macaque**: Falls back to default Ensembl behavior (no override)

---

## Cache Management

All fetched resources are cached:

```bash
~/.cache/sirnaforge/
├── transcriptomes/        # Ensembl cDNA + BWA indices
├── mirna/                 # MirGeneDB/miRBase FASTA files
├── nextflow/              # Pipeline work directories
└── cache.json             # Metadata (TTL, checksums, URLs)
```

### View Cache Status

```bash
sirnaforge cache list
```

### Clear Stale Data

```bash
sirnaforge cache clear --older-than 90
```

---

## Common Workflows

### 1. Human-Only Analysis (Fast)

```bash
sirnaforge workflow TP53 --species human --top-n 20
```

### 2. Multi-Species Comparative

```bash
sirnaforge workflow TP53 --species human,mouse,rat,macaque --top-n 50
```

### 3. Custom Transcriptome + Standard miRNA

```bash
sirnaforge workflow CUSTOM_GENE \
  --input-fasta my_isoforms.fa \
  --transcriptome-fasta my_isoforms.fa \
  --species human
```

**Why repeat the file?**
- `--input-fasta`: Skips gene search, uses these sequences for siRNA design
- `--transcriptome-fasta`: Uses these sequences for off-target checking

### 4. Pre-Indexed Genomes (Production)

```bash
sirnaforge workflow BRCA1 \
  --offtarget-indices \
    human:/mnt/refs/GRCh38,\
    mouse:/mnt/refs/GRCm39,\
    rat:/mnt/refs/rn7
```

---

## Troubleshooting

### "No transcriptome data for species X"

**Cause:** Only 4 species have pre-configured Ensembl sources (human, mouse, rat, macaque)
**Solution:** Provide custom FASTA via `--transcriptome-fasta`

### "Species not recognized" or "Invalid species code"

**Cause:** Rare - species not in registry
**Solution:** Species names are auto-normalized. Supported formats include:
- Common names: `human`, `mouse`, `rat`, `macaque`, `chicken`, `pig`
- miRBase codes: `hsa`, `mmu`, `rno`, `mml`, `gga`, `ssc`
- Scientific names: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`

If you encounter this error, verify spelling or check `src/sirnaforge/data/species_registry.py` for supported species.

### "BWA index not found"

**Cause:** Automatic index building failed or was interrupted
**Solution:** Check cache dir permissions, re-run with `--force-refresh`, or provide `--offtarget-indices`

### Custom FASTA not detected

**Cause:** System treats URLs/paths as single sources
**Solution:** Use comma-separated list (future feature) or concatenate files manually

---

## Implementation Notes

### Current Limitations

1. **Additive transcriptomes**: Not yet exposed via CLI (internal support exists)
2. **Custom miRNA FASTA**: No direct CLI override (use cache manipulation)
3. **Genomic DNA defaults**: None (must provide `--offtarget-indices` explicitly)

### Future Enhancements

- [ ] Comma-separated `--transcriptome-fasta` for additive mode
- [ ] `--custom-mirna` parameter accepting FASTA files
- [ ] Automatic genomic DNA fetching (NCBI/Ensembl FTP)
- [ ] `--add-sequence` for on-the-fly contig injection

---

## Species Normalization

siRNAforge automatically normalizes species names to **canonical form** (`human`, `mouse`, `rat`, `macaque`, `chicken`, `pig`) using the built-in species registry.

### Supported Input Formats

**Common Names** (preferred):
```bash
--species human,mouse,rat,macaque,chicken,pig
```

**miRBase Three-Letter Codes**:
```bash
--species hsa,mmu,rno,mml,gga,ssc
```

**Scientific Names**:
```bash
--species "Homo sapiens,Mus musculus,Rattus norvegicus,Macaca mulatta"
```

**Mixed Formats** (all normalized consistently):
```bash
--species hsa,mouse,"Rattus norvegicus",mml
# Internally stored as: human,mouse,rat,macaque
```

### Why Normalization Matters

Without normalization, you'd have to remember:
- `--species human --mirna-species hsa` (miRNA needs miRBase code)
- Transcriptome expects `"mouse"`, miRNA expects `"mmu"`
- Cross-system comparisons fail due to format mismatch

**With normalization**, all these work identically:
```bash
sirnaforge workflow TP53 --species human
sirnaforge workflow TP53 --species hsa
sirnaforge workflow TP53 --species "Homo sapiens"
```

### Registry Lookup

Species mappings are defined in `src/sirnaforge/data/species_registry.py`:
- **Primary canonical**: `human`, `mouse`, `rat`, `macaque`, `chicken`, `pig`
- **miRBase codes**: `hsa`, `mmu`, `rno`, `mml`, `gga`, `ssc`
- **Aliases**: `rhesus` → `macaque`, `homo_sapiens` → `human`
- **Scientific names**: `Homo sapiens`, `Mus musculus`, `Rattus norvegicus`, etc.

---

## API Usage

For programmatic control, use the Python API:

```python
from sirnaforge.workflow import create_workflow, WorkflowConfig
from sirnaforge.config import ReferenceSelection, ReferenceChoice

# Multi-source transcriptome (internal API)
transcriptome_selection = ReferenceSelection(
    choices=(
        ReferenceChoice.default("ensembl_human_cdna"),
        ReferenceChoice.explicit("/custom/sequences.fa"),
    )
)

workflow = await create_workflow(
    gene_query="TP53",
    output_dir="/results",
    species=["human", "mouse"],
    mirna_database="mirgenedb",
    transcriptome_selection=transcriptome_selection,
)

results = await workflow.run()
```

---

## Summary

| Parameter | Scope | Default | Override Example |
|-----------|-------|---------|------------------|
| `--species` | All layers | 7 species (miRNA) | `--species human,mouse` |
| `--mirna-species` | miRNA only | Maps from `--species` | `--mirna-species hsa,mmu` |
| `--transcriptome-fasta` | Transcriptome only | Ensembl cDNA (4 species) | `--transcriptome-fasta custom.fa` |
| `--offtarget-indices` | BWA indices | Auto-build from transcriptome | `--offtarget-indices human:/data/idx` |

**Key insight:** One parameter (`--species`) intelligently drives everything, with surgical overrides available when needed.
