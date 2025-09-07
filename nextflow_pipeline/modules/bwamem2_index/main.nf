#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Small Nextflow pipeline to build a BWA‑MEM2 index for a given FASTA
 * This is a local helper modeled after nf-core's bwamem2 index module.
 *
 * Usage:
 * nextflow run nextflow_pipeline/modules/bwamem2_index/main.nf --species human --fasta_url "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" --outdir indices
 */

params.species = null
params.fasta_url = null
params.outdir = 'indices'
params.container = ''   // optional docker image to run indexing inside (set to image name to use)

if (!params.species) {
    error "--species is required"
}
if (!params.fasta_url) {
    error "--fasta_url is required (can also be a local path)"
}

process BUILD_BWA_INDEX {
    tag { params.species }
    publishDir "${params.outdir}", mode: 'copy'
    container params.container

    input:
    val species_from_param
    val fasta_url_from_param

    output:
    path "${params.species}_bwa_index", emit: index_dir

    script:
    """
    set -euo pipefail
    mkdir -p ${params.species}_bwa_index
    cd ${params.species}_bwa_index
    # Determine FASTA source: if fasta_url_from_param looks like a URL, download it; otherwise treat it as a local path and copy it
    if echo "${fasta_url_from_param}" | grep -qE '^https?://'; then
        echo "Downloading FASTA from URL: ${fasta_url_from_param}"
        if command -v curl >/dev/null 2>&1; then
            curl -L -o reference.fa.gz "${fasta_url_from_param}"
        else
            wget -O reference.fa.gz "${fasta_url_from_param}"
        fi
        # Decompress if gzipped
        if file reference.fa.gz | grep -q gzip; then
            gunzip -f reference.fa.gz
            mv reference.fa ${params.species}.fa || true
            # If the file had a different name after gunzip, try to normalise
            for f in *.fa *.fasta; do
                if [ -f "\$f" ]; then mv -f "\$f" ${params.species}.fa; break; fi
            done
        else
            mv reference.fa.gz ${params.species}.fa
        fi
    else
        echo "Using local FASTA path: ${fasta_url_from_param}"
        cp "${fasta_url_from_param}" ${params.species}.fa
    fi

    REF=${params.species}.fa
    if [ ! -s "\$REF" ]; then
        echo "ERROR: reference FASTA \$REF not found or empty" >&2
        exit 2
    fi

    echo "Building index for ${params.species} using bwa-mem2 (or bwa fallback)"
    if command -v bwa-mem2 >/dev/null 2>&1; then
        bwa-mem2 index "\$REF"
    elif command -v bwa >/dev/null 2>&1; then
        bwa index "\$REF"
    else
        echo "ERROR: bwa-mem2 or bwa not found in PATH" >&2
        exit 3
    fi

    # Move index files into a directory with prefix name
    mkdir -p ${params.species}_bwa_index
    # BWA creates files with prefix \$REF (filename). Move them into the index dir.
    for f in \$(ls ${params.species}.fa* 2>/dev/null || true); do
        mv -f "\$f" ${params.species}_bwa_index/ || true
    done

    echo "Index built: \$(ls -1 ${params.species}_bwa_index | tr '\\n' ' ')"
    """
}

workflow {
    // If the provided fasta_url points to a local file, pass its absolute path so the process can copy it.
    def fasta_arg = file(params.fasta_url).exists() ? file(params.fasta_url).toAbsolutePath().toString() : params.fasta_url
    BUILD_BWA_INDEX(params.species, fasta_arg)
}
