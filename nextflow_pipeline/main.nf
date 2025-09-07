#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * siRNA Off-target Analysis Pipeline
 * 
 * Clean command-line driven pipeline for siRNA off-target analysis
 * using reusable modules and minimal embedded code.
 */

// Parameters with defaults
params.input = 'nextflow_pipeline/candidates.fasta'                    // Input FASTA with siRNA candidates (default test file)
params.outdir = 'results'             // Output directory
params.genome_species = 'human,rat,rhesus'  // Comma-separated species list
params.sirna_length = 21              // siRNA length
params.max_hits = 10000               // Maximum hits per candidate
params.bwa_k = 12                     // BWA seed length
params.bwa_T = 15                     // BWA minimum score threshold
params.seed_start = 2                 // Seed region start (1-based)
params.seed_end = 8                   // Seed region end (1-based)
params.download_indexes = false      // If true, automatically download s3:// index prefixes using aws CLI
params.genome_indices = [:]           // Map of species -> index paths
params.genome_config = 'nextflow_pipeline/genomes.yaml' // Single genome config file (YAML)

// Initialization is done inside the workflow to satisfy DSL2 requirements

/*
 * Process to validate input sequences using CLI tool
 */
process VALIDATE_INPUT {
    tag "validate"
    publishDir "${params.outdir}/validation", mode: 'copy'
    
    input:
    path input_fasta
    
    output:
    path "validation_report.txt", emit: report
    path "validated_sequences.fasta", emit: sequences
    
    script:
    """
    # Use the siRNA CLI validation command
    uv run sirna validate ${input_fasta} > validation_report.txt 2>&1 || true
    
    # Simple sequence filtering using awk
    awk -v len=${params.sirna_length} '
    BEGIN { valid=0; invalid=0 }
    /^>/ { 
        if (seq != "" && length(seq) == len && seq ~ /^[ATCG]*\$/) {
            print header; print seq; valid++
        } else if (seq != "") {
            invalid++
        }
        header = \$0; seq = ""
    }
    !/^>/ { seq = seq \$0 }
    END { 
        if (seq != "" && length(seq) == len && seq ~ /^[ATCG]*\$/) {
            print header; print seq; valid++
        }
        print "Validated " valid " sequences, rejected " invalid > "/dev/stderr"
    }' ${input_fasta} > validated_sequences.fasta
    
    echo "Sequence validation completed" >> validation_report.txt
    echo "Valid sequences: \$(grep -c '^>' validated_sequences.fasta)" >> validation_report.txt
    # Fail if no valid sequences found
    if [ \$(grep -c '^>' validated_sequences.fasta) -eq 0 ]; then
        echo "ERROR: No valid sequences found after validation" >&2
        exit 1
    fi
    """
}

/*
 * Process to run off-target analysis using the wrapper CLI
 */
process RUN_OFFTARGET_ANALYSIS {
    tag "${species}"
    publishDir "${params.outdir}/${species}", mode: 'copy'
    
    input:
    val species
    path sequences
    
    output:
    path "${species}_offtargets.tsv", emit: tsv
    path "${species}_offtargets.json", emit: json
    path "${species}_summary.txt", emit: summary
    
    script:
    // Determine index path and optional container from config at runtime (resolve inside shell)
    def base_index = "/data/genomes/${species}/index"
    def container_image = null

    // Build the wrapper command using a runtime INDEX_PATH environment variable so it can be swapped after download
    def base_cmd = "python3 ${moduleDir}/offtarget_wrapper.py --mode combined --index-prefix \"\$INDEX_PATH\" --queries ${sequences} --out-prefix ${species}_offtargets --bwa-k ${params.bwa_k} --bwa-T ${params.bwa_T} --max-hits ${params.max_hits} --seed-start ${params.seed_start} --seed-end ${params.seed_end} --verbose"
    // If running inside a container, include mounts; escape $ so it is evaluated at runtime inside the container
    def wrapper_cmd = container_image ? "docker run --rm -v \$PWD:/work -v ${moduleDir}:${moduleDir} -w /work ${container_image} /bin/bash -c \"${base_cmd}\"" : base_cmd

    """
    # Resolve index path using Python (reads genomes.yaml if present, else uses --genome_indices)
    INDEX_PATH="${base_index}"
    IDX=\$(python3 ${moduleDir}/scripts/resolve_index.py "${species}" "${params.genome_config}")
    if [ -n "\$IDX" ]; then
        INDEX_PATH="\$IDX"
    fi

    # If index is an s3:// path and downloads are enabled, try to fetch it locally
    if [[ "\${INDEX_PATH}" == s3://* ]]; then
        if [ "${params.download_indexes}" = "true" ] || [ "\${DOWNLOAD_INDEXES}" = "true" ]; then
            if command -v aws >/dev/null 2>&1; then
                TMP_INDEX_DIR=".downloaded_indexes/${species}"
                mkdir -p "\${TMP_INDEX_DIR}"
                echo "Downloading \${INDEX_PATH} -> \${TMP_INDEX_DIR}" >&2
                aws s3 --no-sign-request sync "\${INDEX_PATH}" "\${TMP_INDEX_DIR}" || echo "Warning: aws s3 sync failed for \${INDEX_PATH}" >&2
                # prefer prefix inside tmp dir
                INDEX_PATH="\${TMP_INDEX_DIR}"

                    # If downloaded content looks like a reference FASTA, build a bwa index
                    # Look for common FASTA names (.fa, .fasta, .fa.gz)
                    FASTA_FILE=""
                    # prefer uncompressed fasta first
                    for f in "\${TMP_INDEX_DIR}"/*.fa "\${TMP_INDEX_DIR}"/*.fasta; do
                        if [ -f "\$f" ]; then FASTA_FILE="\$f"; break; fi
                    done
                    # then gzipped fasta
                    if [ -z "\$FASTA_FILE" ]; then
                        for f in "\${TMP_INDEX_DIR}"/*.fa.gz "\${TMP_INDEX_DIR}"/*.fasta.gz; do
                            if [ -f "\$f" ]; then FASTA_FILE="\$f"; break; fi
                        done
                        # decompress if needed
                        if [ -n "\$FASTA_FILE" ] && file "\$FASTA_FILE" | grep -q gzip; then
                            echo "Decompressing \$FASTA_FILE" >&2
                            gunzip -f "\$FASTA_FILE"
                            # update FASTA_FILE to decompressed name
                            FASTA_FILE="\${FASTA_FILE%.gz}"
                        fi
                    fi

                    if [ -n "\$FASTA_FILE" ]; then
                        echo "Found reference FASTA: \$FASTA_FILE. Building index in \${TMP_INDEX_DIR}/index" >&2
                        # prefer helper script bundled in the pipeline
                        if [ -x "\${moduleDir}/scripts/build_index.sh" ]; then
                            bash "\${moduleDir}/scripts/build_index.sh" "\$FASTA_FILE" "\${TMP_INDEX_DIR}/index" || echo "Warning: index build failed" >&2
                            INDEX_PATH="\${TMP_INDEX_DIR}/index"
                        else
                            # fallback to bwa-mem2 or bwa if available
                            if command -v bwa-mem2 >/dev/null 2>&1; then
                                bwa-mem2 index -p "\${TMP_INDEX_DIR}/index" "\$FASTA_FILE" || echo "Warning: bwa-mem2 index failed" >&2
                                INDEX_PATH="\${TMP_INDEX_DIR}/index"
                            elif command -v bwa >/dev/null 2>&1; then
                                bwa index -p "\${TMP_INDEX_DIR}/index" "\$FASTA_FILE" || echo "Warning: bwa index failed" >&2
                                INDEX_PATH="\${TMP_INDEX_DIR}/index"
                            else
                                echo "Warning: no index builder (build_index.sh, bwa-mem2 or bwa) found; leaving INDEX_PATH=\${TMP_INDEX_DIR}" >&2
                            fi
                        fi
                    fi
            else
                echo "Warning: aws CLI not found; cannot download \${INDEX_PATH}" >&2
            fi
        else
            echo "Note: index \${INDEX_PATH} is an s3:// path; set --download_indexes true to auto-download" >&2
        fi
    fi

    # Check if index exists (local or downloaded)
    if [ ! -f "\${INDEX_PATH}.fasta" ] && [ ! -f "\${INDEX_PATH}.fa" ] && [ ! -d "\${INDEX_PATH}" ]; then
        echo "Warning: Index \${INDEX_PATH} not found for ${species}" >&2
        
        # Create empty results
        echo -e "qname\\tqseq\\trname\\tcoord\\tstrand\\tcigar\\tmapq\\tAS\\tNM\\tmismatch_count\\tmismatch_positions\\tofftarget_score\\tsource" > ${species}_offtargets.tsv
        echo "[]" > ${species}_offtargets.json
        
        echo "Off-target Analysis Summary for ${species}" > ${species}_summary.txt
        echo "Status: Index not found" >> ${species}_summary.txt
        echo "Total hits: 0" >> ${species}_summary.txt
        
        exit 0
    fi
    
    # If the validated sequences file is empty, skip this species
    if [ ! -f "${sequences}" ] || [ \$(grep -c '^>' ${sequences}) -eq 0 ]; then
        echo "Warning: No validated sequences to analyze for ${species}" >&2
        echo -e "qname\tqseq\trname\tcoord\tstrand\tcigar\tmapq\tAS\tNM\tmismatch_count\tmismatch_positions\tofftarget_score\tsource" > ${species}_offtargets.tsv
        echo "[]" > ${species}_offtargets.json
        echo "Off-target Analysis Summary for ${species}" > ${species}_summary.txt
        echo "Status: No validated sequences" >> ${species}_summary.txt
        echo "Total hits: 0" >> ${species}_summary.txt
        exit 0
    fi
    # Export INDEX_PATH for wrapper command and run off-target analysis using the Python wrapper (or in container)
    export INDEX_PATH
    ${wrapper_cmd} || {
        
        # Handle failure gracefully
        echo "Error: Off-target analysis failed for ${species}" >&2
        
        # Create empty results on failure
        echo -e "qname\\tqseq\\trname\\tcoord\\tstrand\\tcigar\\tmapq\\tAS\\tNM\\tmismatch_count\\tmismatch_positions\\tofftarget_score\\tsource" > ${species}_offtargets.tsv
        echo "[]" > ${species}_offtargets.json
        
        echo "Off-target Analysis Summary for ${species}" > ${species}_summary.txt
        echo "Status: Analysis failed" >> ${species}_summary.txt
        echo "Total hits: 0" >> ${species}_summary.txt
        
        exit 0
    }
    
    # Generate summary using simple shell commands
    TOTAL_HITS=\$(tail -n +2 ${species}_offtargets.tsv | wc -l)
    
    echo "Off-target Analysis Summary for ${species}" > ${species}_summary.txt
    echo "Status: Completed successfully" >> ${species}_summary.txt
    echo "Total hits: \$TOTAL_HITS" >> ${species}_summary.txt
    echo "Index used: \${INDEX_PATH}" >> ${species}_summary.txt
    echo "Parameters: k=${params.bwa_k}, T=${params.bwa_T}" >> ${species}_summary.txt
    """
}

/*
 * Process to aggregate results using shell commands and simple tools
 */
process AGGREGATE_RESULTS {
    tag "aggregate"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path tsv_files
    path json_files
    path summary_files
    
    output:
    path "combined_offtargets.tsv", emit: combined_tsv
    path "combined_offtargets.json", emit: combined_json
    path "final_summary.txt", emit: final_summary
    path "offtarget_report.html", emit: html_report
    
    script:
    """
    # Combine TSV files using standard Unix tools
    echo "Combining TSV results..."
    
    # Get header from first file
    head -n1 \$(ls *.tsv | head -n1) > combined_offtargets.tsv
    
    # Append data from all files, adding species column
    for tsv_file in *.tsv; do
        species=\$(basename "\$tsv_file" "_offtargets.tsv")
        tail -n +2 "\$tsv_file" | awk -v sp="\$species" '{print \$0 "\\t" sp}' >> combined_offtargets.tsv
    done
    
    # Fix header to include species column
    sed -i '1s/\$/\\tspecies/' combined_offtargets.tsv
    
    # Combine JSON files using jq if available, otherwise use Python one-liner
    echo "Combining JSON results..."
    
    if command -v jq >/dev/null 2>&1; then
        # Use jq for JSON manipulation
        jq -s 'add' *.json > combined_offtargets.json
    else
        # Fallback to Python one-liner
        python3 -c "
import json, glob, sys
combined = []
for f in glob.glob('*_offtargets.json'):
    species = f.replace('_offtargets.json', '')
    try:
        with open(f) as fh:
            data = json.load(fh)
            for item in data:
                item['species'] = species
            combined.extend(data)
    except: pass
with open('combined_offtargets.json', 'w') as out:
    json.dump(combined, out, indent=2)
"
    fi
    
    # Generate summary using shell commands
    echo "Generating final summary..."
    
    TOTAL_HITS=\$(tail -n +2 combined_offtargets.tsv | wc -l)
    UNIQUE_QUERIES=\$(tail -n +2 combined_offtargets.tsv | cut -f1 | sort -u | wc -l)
    SPECIES_COUNT=\$(tail -n +2 combined_offtargets.tsv | cut -f13 | sort -u | wc -l)
    
    cat > final_summary.txt << EOF
siRNA Off-target Analysis - Final Summary
========================================

Analysis Parameters:
  Input file: ${params.input}
  Genome species: ${params.genome_species}
  siRNA length: ${params.sirna_length}
  Max hits per query: ${params.max_hits}

Overall Results:
  Total hits across all species: \$TOTAL_HITS
  Unique query sequences: \$UNIQUE_QUERIES
  Species analyzed: \$SPECIES_COUNT

Results by Species:
EOF

    # Add per-species counts
    for summary_file in *_summary.txt; do
        if [[ "\$summary_file" != "final_summary.txt" ]]; then
            species=\$(basename "\$summary_file" "_summary.txt")
            hits=\$(grep "Total hits:" "\$summary_file" | cut -d: -f2 | tr -d ' ')
            echo "  \$species: \$hits hits" >> final_summary.txt
        fi
    done
    
    # Generate simple HTML report
    cat > offtarget_report.html << EOF
<!DOCTYPE html>
<html>
<head>
    <title>siRNA Off-target Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        .summary { background-color: #f9f9f9; padding: 15px; border-left: 5px solid #4CAF50; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 siRNA Off-target Analysis Report</h1>
        <p>Generated on: \$(date)</p>
    </div>
    
    <div class="section">
        <h2>Analysis Parameters</h2>
        <ul>
            <li><strong>Input file:</strong> ${params.input}</li>
            <li><strong>Genome species:</strong> ${params.genome_species}</li>
            <li><strong>siRNA length:</strong> ${params.sirna_length} nt</li>
            <li><strong>Max hits per query:</strong> ${params.max_hits}</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Results Summary</h2>
        <div class="summary">
            <p><strong>Total off-target hits:</strong> \$TOTAL_HITS</p>
            <p><strong>Unique queries:</strong> \$UNIQUE_QUERIES</p>
            <p><strong>Species analyzed:</strong> \$SPECIES_COUNT</p>
        </div>
        
        <h3>Hits by Species</h3>
        <table>
            <tr><th>Species</th><th>Number of Hits</th></tr>
EOF

    # Add species data to HTML table
    for summary_file in *_summary.txt; do
        if [[ "\$summary_file" != "final_summary.txt" ]]; then
            species=\$(basename "\$summary_file" "_summary.txt")
            hits=\$(grep "Total hits:" "\$summary_file" | cut -d: -f2 | tr -d ' ')
            echo "            <tr><td>\$species</td><td>\$hits</td></tr>" >> offtarget_report.html
        fi
    done
    
    cat >> offtarget_report.html << 'EOF'
        </table>
    </div>
    
    <div class="section">
        <h2>Output Files</h2>
        <ul>
            <li><strong>combined_offtargets.tsv:</strong> All hits in tabular format</li>
            <li><strong>combined_offtargets.json:</strong> All hits in JSON format</li>
            <li><strong>final_summary.txt:</strong> Text summary of results</li>
            <li><strong>[species]_offtargets.*:</strong> Species-specific results</li>
        </ul>
    </div>
</body>
</html>
EOF
    
    echo "Results aggregation completed successfully"
    """
}

/*
 * Main workflow
 */
workflow {
    // Parameter validation and logging
    if (!params.input) {
        error "Input FASTA file is required. Use --input <file>"
    }

    log.info """\
    siRNA OFF-TARGET ANALYSIS PIPELINE
    ===================================
    Input FASTA     : ${params.input}
    Output directory: ${params.outdir}
    Genome species  : ${params.genome_species}
    siRNA length    : ${params.sirna_length}
    Max hits        : ${params.max_hits}
    BWA seed length : ${params.bwa_k}
    BWA threshold   : ${params.bwa_T}
    Seed region     : ${params.seed_start}-${params.seed_end}
    """

    // genomes parsed at top-level and exposed via params.genomes

    // Parse species list
    def species_list = params.genome_species.split(',').collect { it.trim() }
    species_ch = Channel.fromList(species_list)
    
    // Input validation
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Validate input sequences
    VALIDATE_INPUT(input_ch)

    // Run off-target analysis for each species (genomes map is available via params.genomes)
    RUN_OFFTARGET_ANALYSIS(
        species_ch,
        VALIDATE_INPUT.out.sequences
    )
    
    // Aggregate results
    AGGREGATE_RESULTS(
        RUN_OFFTARGET_ANALYSIS.out.tsv.collect(),
        RUN_OFFTARGET_ANALYSIS.out.json.collect(),
        RUN_OFFTARGET_ANALYSIS.out.summary.collect()
    )
}

/*
 * Workflow completion message handled by final aggregation and Nextflow reports
 */
