# ZFN reference fixtures

`ccr5_grch38_on_target.fa` is an 89 bp slice of the GRCh38 primary-assembly
chromosome 3 record, 1-based inclusive coordinates 46,373,023–46,373,111.
Its 29 bp center is the CCR5 on-target sequence published by Fine et al.:
`GTCATCCTCATCCTGATAAACTGCAAAAG`. The exact sequence is present in GRCh38
chr3:46,373,053–46,373,081. The current GRCh38 chromosome record is
[NC_000003.12](https://www.ncbi.nlm.nih.gov/nuccore/NC_000003.12).

The source FASTA was the Ensembl-style record
`>3 dna:chromosome chromosome:GRCh38:3:1:198295559:1 REF`; its SHA-256 was
`1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b`.
Only this tiny locus slice is checked in. This fixture verifies sequence and
query orientation on the reference; the PROGNOS table's coordinate 46,414,544
is hg19 and must not be reported as a GRCh38 coordinate.

Source sequence: Fine et al. (2014), *An online bioinformatics tool predicts
zinc finger and TALE nuclease off-target cleavage*,
[doi:10.1093/nar/gkt1326](https://doi.org/10.1093/nar/gkt1326).
