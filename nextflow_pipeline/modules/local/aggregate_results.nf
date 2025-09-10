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
from sirnaforge.core.off_target import parse_fasta_file, write_fasta_file
from pathlib import Path
import glob
import json

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print('Warning: pandas not available, using basic aggregation')

# Collect all analysis and summary files
analysis_files = glob.glob('*.tsv')
summary_files = glob.glob('*_summary.txt')

print(f'Found {len(analysis_files)} analysis files and {len(summary_files)} summary files')

# Combine all TSV analysis results
if analysis_files and PANDAS_AVAILABLE:
    all_results = []
    for tsv_file in analysis_files:
        try:
            df = pd.read_csv(tsv_file, sep='\\t')
            df['source_file'] = tsv_file
            all_results.append(df)
        except Exception as e:
            print(f'Warning: Could not read {tsv_file}: {e}')

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df.to_csv('combined_offtarget_analysis.tsv', sep='\\t', index=False)
        print(f'Combined {len(all_results)} analysis files into combined_offtarget_analysis.tsv')
elif analysis_files:
    # Basic aggregation without pandas
    with open('combined_offtarget_analysis.tsv', 'w') as outf:
        header_written = False
        for tsv_file in analysis_files:
            try:
                with open(tsv_file, 'r') as inf:
                    lines = inf.readlines()
                    if lines:
                        if not header_written:
                            outf.write(lines[0])  # Write header
                            header_written = True
                        for line in lines[1:]:  # Skip header in subsequent files
                            outf.write(line)
            except Exception as e:
                print(f'Warning: Could not read {tsv_file}: {e}')

# Aggregate summary information
import datetime
summary_data = {
    'genome_species': '${genome_species}'.split(','),
    'total_analysis_files': len(analysis_files),
    'total_summary_files': len(summary_files),
    'analysis_timestamp': str(datetime.datetime.now()),
    'file_details': []
}

for summary_file in summary_files:
    try:
        with open(summary_file, 'r') as f:
            content = f.read()
            summary_data['file_details'].append({
                'file': summary_file,
                'content': content
            })
    except Exception as e:
        print(f'Warning: Could not read {summary_file}: {e}')

# Write combined summary
with open('combined_summary.json', 'w') as f:
    json.dump(summary_data, f, indent=2)

# Write final summary
with open('final_summary.txt', 'w') as f:
    f.write('SiRNA Off-Target Analysis Summary\\n')
    f.write('=' * 40 + '\\n')
    f.write(f'Genome species analyzed: {len(summary_data[\"genome_species\"])}\\n')
    f.write(f'Analysis files processed: {len(analysis_files)}\\n')
    f.write(f'Summary files processed: {len(summary_files)}\\n')
    f.write(f'Analysis completed at: {summary_data[\"analysis_timestamp\"]}\\n')

    if analysis_files:
        total_hits = 0
        try:
            with open('combined_offtarget_analysis.tsv', 'r') as f:
                total_hits = len(f.readlines()) - 1  # Subtract header
            f.write(f'Total off-target hits found: {total_hits}\\n')
        except:
            f.write('Could not calculate total hits\\n')

print('Aggregation completed successfully')
"

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
