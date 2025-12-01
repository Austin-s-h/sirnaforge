# CLI Reference

> **Auto-generated**: All command output below is captured live during documentation build from the actual `sirnaforge` CLI.

This reference shows each command with its real `--help` output and working examples.

## Help & Version

### Main Help

```{program-output} uv run sirnaforge --help
```

### Version

```{program-output} uv run sirnaforge version
```

---

## workflow

Run complete siRNA design from gene query to scored candidates.

### Help

```{program-output} uv run sirnaforge workflow --help
```

:::{note}
The workflow command searches for gene transcripts, designs siRNA candidates, scores them using thermodynamic analysis, and outputs ranked results.
:::

---

## search

Search gene databases and retrieve transcript sequences.

### Help

```{program-output} uv run sirnaforge search --help
```

---

## design

Design siRNA candidates from FASTA sequences.

### Help

```{program-output} uv run sirnaforge design --help
```

### Example: Design from Sample Data

```{program-output} uv run sirnaforge design ../examples/sample_transcripts.fasta -o /tmp/sirna_example.csv --top-n 5
```

#### Output Preview

```{program-output} head -6 /tmp/sirna_example.csv
:shell:
```

---

## validate

Check FASTA file format and content.

### Help

```{program-output} uv run sirnaforge validate --help
```

### Example: Validate Sample Data

```{program-output} uv run sirnaforge validate ../examples/sample_transcripts.fasta
```

---

## config

Show default configuration parameters.

```{program-output} uv run sirnaforge config
```

---

## sequences

Manage siRNA sequences and chemical modification metadata.

### Help

```{program-output} uv run sirnaforge sequences --help
```

---

## cache

Manage miRNA database cache for off-target analysis.

### Help

```{program-output} uv run sirnaforge cache --help
```
