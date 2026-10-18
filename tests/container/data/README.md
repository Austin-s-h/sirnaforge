# Chromosome-scoped CI reference

`human_grch38_chr21.cdna.fa.gz` contains **all 4,564 chromosome 21 cDNA records**
from a cached Ensembl human GRCh38 snapshot, preserving headers, sequences and
source order. It has 12,534,162 sequence bases and compresses to about 1.5 MB.
The adjacent JSON records the source URL, download timestamp, source SHA-256,
selection rule and fixture checksums. The cache did not record a release number;
the mutable `current_fasta` URL must not be treated as a release identifier.

The container test decompresses and verifies this fixture, passes it through
`--transcriptome-fasta`, builds a fresh BWA-MEM2 index and requires completed
transcriptome evidence and a rendered report. This replaces a full mouse cDNA
download/index build that was killed on the GitHub runner. It validates packaging
and the custom-reference path within CI resources; screening counts are scoped
to chromosome 21.
Both transcriptome reference tests seed their human miRNA cache from the existing
toy fixture, so their success does not depend on a live miRNA service either.

To regenerate from the exact **uncompressed** source identified in the JSON:

```python
from pathlib import Path
import gzip
import hashlib
import json

fixture = Path("tests/container/data/human_grch38_chr21.cdna.fa.gz")
metadata = json.loads(fixture.with_suffix(".json").read_text())
source = Path("Homo_sapiens.GRCh38.cdna.all.fa")
digest = hashlib.sha256()
selected = bytearray()
keep = False
with source.open("rb") as handle:
    for line in handle:
        digest.update(line)
        if line.startswith(b">"):
            keep = b" chromosome:GRCh38:21:" in line
        if keep:
            selected.extend(line)
assert digest.hexdigest() == metadata["source_uncompressed_sha256"]
assert hashlib.sha256(selected).hexdigest() == metadata["fasta_sha256"]
fixture.write_bytes(gzip.compress(bytes(selected), mtime=0))
```

The uncompressed checksum is independent of Python/zlib gzip implementation
details. Deliberately changing the source snapshot requires updating the JSON.
