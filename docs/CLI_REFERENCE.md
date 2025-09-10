# 🧬 siRNAforge CLI Reference

```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge [OPTIONS] COMMAND [ARGS]...[0m[1m                                                                         [0m[1m [0m
[1m                                                                                                                        [0m
 🧬 siRNAforge - Comprehensive siRNA design toolkit for gene silencing                                                  
                                                                                                                        
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-install[0m[1;36m-completion[0m          Install completion for the current shell.                                              [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-show[0m[1;36m-completion[0m             Show completion for the current shell, to copy it or customize the installation.       [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m                        Show this message and exit.                                                            [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m
[2m╭─[0m[2m Commands [0m[2m──────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36msearch   [0m[1;36m [0m 🔍 Search for gene transcripts and retrieve sequences.                                                    [2m│[0m
[2m│[0m [1;36mworkflow [0m[1;36m [0m 🧬 Run complete siRNA design workflow from gene query to off-target analysis.                             [2m│[0m
[2m│[0m [1;36mdesign   [0m[1;36m [0m 🎯 Design siRNA candidates from transcript sequences.                                                     [2m│[0m
[2m│[0m [1;36mvalidate [0m[1;36m [0m 🔍 Validate input FASTA file format and content.                                                          [2m│[0m
[2m│[0m [1;36mversion  [0m[1;36m [0m 📦 Show version information.                                                                              [2m│[0m
[2m│[0m [1;36mconfig   [0m[1;36m [0m ⚙️  Show default configuration parameters.                                                                 [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `search`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge search [OPTIONS] QUERY[0m[1m                                                                              [0m[1m [0m
[1m                                                                                                                        [0m
 🔍 Search for gene transcripts and retrieve sequences.                                                                 
                                                                                                                        
[2m╭─[0m[2m Arguments [0m[2m─────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [31m*[0m    query      [1;33mTEXT[0m  Gene ID, gene name, or transcript ID to search for [2;31m[required][0m                                  [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-output[0m             [1;32m-o[0m                            [1;33mPATH[0m  Output FASTA file for transcript sequences                  [2m│[0m
[2m│[0m                                                          [2m[default: transcripts.fasta]              [0m                  [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-database[0m           [1;32m-d[0m                            [1;33mTEXT[0m  Database to search (ensembl, refseq, gencode)               [2m│[0m
[2m│[0m                                                          [2m[default: ensembl]                           [0m               [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-all[0m                [1;32m-a[0m                            [1;33m    [0m  Search all databases                                        [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-no[0m[1;36m-sequence[0m                                      [1;33m    [0m  Skip sequence retrieval (metadata only)                     [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-canonical[0m[1;36m-only[0m                                   [1;33m    [0m  Extract only canonical isoforms                             [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-extract[0m[1;36m-canonical[0m      [1;35m-[0m[1;35m-no[0m[1;35m-extract-canonical[0m    [1;33m    [0m  Automatically extract canonical isoforms to separate file   [2m│[0m
[2m│[0m                                                          [2m[default: extract-canonical]                             [0m   [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-types[0m              [1;32m-t[0m                            [1;33mTEXT[0m  Comma-separated list of transcript types to include (e.g.,  [2m│[0m
[2m│[0m                                                          protein_coding,lncRNA)                                      [2m│[0m
[2m│[0m                                                          [2m[default: protein_coding,lncRNA]                           [0m [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-exclude[0m[1;36m-types[0m                                    [1;33mTEXT[0m  Comma-separated list of transcript types to exclude         [2m│[0m
[2m│[0m                                                          [2m[default: nonsense_mediated_decay,retained_intron] [0m         [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-verbose[0m            [1;32m-v[0m                            [1;33m    [0m  Enable verbose output                                       [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m                                             [1;33m    [0m  Show this message and exit.                                 [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `workflow`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge workflow [OPTIONS] GENE_QUERY[0m[1m                                                                       [0m[1m [0m
[1m                                                                                                                        [0m
 🧬 Run complete siRNA design workflow from gene query to off-target analysis.                                          
                                                                                                                        
[2m╭─[0m[2m Arguments [0m[2m─────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [31m*[0m    gene_query      [1;33mTEXT[0m  Gene name or ID to analyze [2;31m[required][0m                                                     [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-input[0m[1;36m-fasta[0m             [1;33mFILE                       [0m  Path to an input FASTA file to use instead of performing a    [2m│[0m
[2m│[0m                                                        gene search                                                   [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-output[0m[1;36m-dir[0m      [1;32m-o[0m      [1;33mPATH                       [0m  Output directory for all workflow results                     [2m│[0m
[2m│[0m                                                        [2m[default: sirna_workflow_output]         [0m                     [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-database[0m        [1;32m-d[0m      [1;33mTEXT                       [0m  Database to search (ensembl, refseq, gencode)                 [2m│[0m
[2m│[0m                                                        [2m[default: ensembl]                           [0m                 [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-top[0m[1;36m-n[0m           [1;32m-n[0m      [1;33mINTEGER RANGE [x>=1[0m[1;2;33m][0m[1;33m       [0m  Number of top siRNA candidates to generate [2m[default: 20][0m      [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-offtarget[0m[1;36m-n[0m             [1;33mINTEGER RANGE [x>=1[0m[1;2;33m][0m[1;33m       [0m  Number of top candidates for off-target analysis              [2m│[0m
[2m│[0m                                                        [2m[default: 10]                                   [0m              [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-genome[0m[1;36m-species[0m          [1;33mTEXT                       [0m  Comma-separated list of genome species for off-target         [2m│[0m
[2m│[0m                                                        analysis                                                      [2m│[0m
[2m│[0m                                                        [2m[default: human,rat,rhesus]                                  [0m [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-gc[0m[1;36m-min[0m                  [1;33mFLOAT RANGE [0.0<=x<=100.0[0m[1;2;33m][0m  Minimum GC content percentage [2m[default: 30.0][0m                 [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-gc[0m[1;36m-max[0m                  [1;33mFLOAT RANGE [0.0<=x<=100.0[0m[1;2;33m][0m  Maximum GC content percentage [2m[default: 52.0][0m                 [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-length[0m          [1;32m-l[0m      [1;33mINTEGER RANGE [19<=x<=23[0m[1;2;33m][0m[1;33m  [0m  siRNA length in nucleotides [2m[default: 21][0m                     [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-verbose[0m         [1;32m-v[0m      [1;33m                           [0m  Enable verbose output                                         [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-log[0m[1;36m-file[0m                [1;33mPATH                       [0m  Path to centralized log file (overrides SIRNAFORGE_LOG_FILE   [2m│[0m
[2m│[0m                                                        env)                                                          [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m                    [1;33m                           [0m  Show this message and exit.                                   [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `design`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge design [OPTIONS] INPUT_FILE[0m[1m                                                                         [0m[1m [0m
[1m                                                                                                                        [0m
 🎯 Design siRNA candidates from transcript sequences.                                                                  
                                                                                                                        
[2m╭─[0m[2m Arguments [0m[2m─────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [31m*[0m    input_file      [1;33mFILE[0m  Input FASTA file containing transcript sequences [2;31m[required][0m                               [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-output[0m            [1;32m-o[0m      [1;33mPATH                       [0m  Output file for siRNA candidates                            [2m│[0m
[2m│[0m                                                          [2m[default: sirna_results.tsv]    [0m                            [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-length[0m            [1;32m-l[0m      [1;33mINTEGER RANGE [19<=x<=23[0m[1;2;33m][0m[1;33m  [0m  siRNA length in nucleotides [2m[default: 21][0m                   [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-top[0m[1;36m-n[0m             [1;32m-n[0m      [1;33mINTEGER RANGE [x>=1[0m[1;2;33m][0m[1;33m       [0m  Number of top candidates to return [2m[default: 10][0m            [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-gc[0m[1;36m-min[0m                    [1;33mFLOAT RANGE [0.0<=x<=100.0[0m[1;2;33m][0m  Minimum GC content percentage [2m[default: 30.0][0m               [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-gc[0m[1;36m-max[0m                    [1;33mFLOAT RANGE [0.0<=x<=100.0[0m[1;2;33m][0m  Maximum GC content percentage [2m[default: 52.0][0m               [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-max[0m[1;36m-poly-runs[0m             [1;33mINTEGER RANGE [x>=1[0m[1;2;33m][0m[1;33m       [0m  Maximum consecutive identical nucleotides [2m[default: 3][0m      [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-genome[0m[1;36m-index[0m              [1;33mPATH                       [0m  Genome index for off-target analysis                        [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-snp[0m[1;36m-file[0m                  [1;33mPATH                       [0m  VCF file with SNPs to avoid                                 [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-skip[0m[1;36m-structure[0m            [1;33m                           [0m  Skip secondary structure prediction (faster)                [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-skip[0m[1;36m-off-targets[0m          [1;33m                           [0m  Skip off-target analysis (faster)                           [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-verbose[0m           [1;32m-v[0m      [1;33m                           [0m  Enable verbose output                                       [2m│[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m                      [1;33m                           [0m  Show this message and exit.                                 [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `validate`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge validate [OPTIONS] INPUT_FILE[0m[1m                                                                       [0m[1m [0m
[1m                                                                                                                        [0m
 🔍 Validate input FASTA file format and content.                                                                       
                                                                                                                        
[2m╭─[0m[2m Arguments [0m[2m─────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [31m*[0m    input_file      [1;33mFILE[0m  FASTA file to validate [2;31m[required][0m                                                         [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m          Show this message and exit.                                                                          [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `config`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge config [OPTIONS][0m[1m                                                                                    [0m[1m [0m
[1m                                                                                                                        [0m
 ⚙️  Show default configuration parameters.                                                                              
                                                                                                                        
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m          Show this message and exit.                                                                          [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```

### `version`
```bash
[1m                                                                                                                        [0m
[1m [0m[1;33mUsage: [0m[1msirnaforge version [OPTIONS][0m[1m                                                                                   [0m[1m [0m
[1m                                                                                                                        [0m
 📦 Show version information.                                                                                           
                                                                                                                        
[2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────────────────────────────────────────────[0m[2m─╮[0m
[2m│[0m [1;36m-[0m[1;36m-help[0m          Show this message and exit.                                                                          [2m│[0m
[2m╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯[0m

```
