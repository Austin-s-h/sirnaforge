#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root inside the local sirnaforge benchmark image.
# The output root is deliberately separate from the historical 4.0.0 workflow results.
export PYTHONPATH="/workspace/src${PYTHONPATH:+:${PYTHONPATH}}"

output_root="workflow_results/oligogym_native_design_preproduction_v1"
reference_root="tests/data/benchmarks/oligogym_native_design/references/representative"
transcriptomes="human:/workspace/${reference_root}/human_protein_coding_gene_representatives.fa,mouse:/workspace/${reference_root}/mouse_protein_coding_gene_representatives.fa,rat:/workspace/${reference_root}/rat_protein_coding_gene_representatives.fa,macaque:/workspace/${reference_root}/macaque_protein_coding_gene_representatives.fa"

run_panel() {
    local dataset="$1"
    local length="$2"
    local input="tests/data/benchmarks/oligogym_native_design/inputs/${dataset}_native_len${length}.fasta"
    local output="${output_root}/${dataset}_len${length}"

    echo "Running ${dataset} sequentially with the pre-production scoring defaults"
    sirnaforge workflow "benchmark_oligogym_${dataset}" \
        --input-fasta "/workspace/${input}" \
        --length "${length}" \
        --run-mode qualified \
        --design-mode sirna \
        --database ensembl \
        --species human,mouse,rat,macaque \
        --query-species human \
        --ortholog-mapping "/workspace/tests/data/benchmarks/oligogym_native_design/orthology/${dataset}.json" \
        --transcriptome-fasta "${transcriptomes}" \
        --gc-min 0 \
        --gc-max 100 \
        --modifications none \
        --output-dir "/workspace/${output}"
}

run_panel ichihara_2007_1 19
run_panel ichihara_2007_2 19
run_panel martinelli_2023_1 21

echo "Completed all three pre-production OligoGym panels under ${output_root}"
