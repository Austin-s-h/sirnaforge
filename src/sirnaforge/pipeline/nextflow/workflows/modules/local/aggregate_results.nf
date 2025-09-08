process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    path analysis_files
    path summary_files
    val genome_species

    output:
    path "combined_*.tsv", emit: combined_analyses, optional: true
    path "combined_summary.json", emit: combined_summary, optional: true
    path "final_summary.txt", emit: final_summary
    path "analysis_report.html", emit: html_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import aggregate_offtarget_results
from pathlib import Path
import glob

# Collect all analysis and summary files
analysis_files = glob.glob('*_analysis.tsv')
summary_files = glob.glob('*_summary.json')

print(f'Found {len(analysis_files)} analysis files and {len(summary_files)} summary files')

if analysis_files or summary_files:
    # Create a temporary results directory structure
    results_dir = Path('temp_results')
    results_dir.mkdir(exist_ok=True)

    # Organize files by species (extract from filename)
    species_list = '${genome_species}'.split(',')
    for species in species_list:
        species_dir = results_dir / species.strip()
        species_dir.mkdir(exist_ok=True)

        # Move relevant files to species directory
        for f in analysis_files:
            if species.strip() in f:
                import shutil
                shutil.copy(f, species_dir / f)

        for f in summary_files:
            if species.strip() in f:
                import shutil
                shutil.copy(f, species_dir / f)

    # Run aggregation
    output_dir = aggregate_offtarget_results(
        results_dir=str(results_dir),
        output_dir='.',
        genome_species='${genome_species}'
    )

    print(f'Aggregation completed: {output_dir}')
else:
    # Create empty final summary
    with open('final_summary.txt', 'w') as f:
        f.write('No analysis results found to aggregate\\n')
    print('No files to aggregate')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch combined_mirna_analysis.tsv
    touch combined_transcriptome_analysis.tsv
    echo '{}' > combined_summary.json
    echo 'Aggregation completed' > final_summary.txt
    touch analysis_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
