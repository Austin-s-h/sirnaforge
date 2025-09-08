process RESOLVE_INDICES {
    tag "$genome_species"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:4ec8e5b31fe18bb4':
        'community.wave.seqera.io/library/python_pyyaml:4ec8e5b31fe18bb4' }"

    input:
    val genome_species
    path genome_config
    val genome_indices_override
    val download_indexes

    output:
    path "resolved_indices.json", emit: indices
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import resolve_genome_indices
import json

# Resolve genome indices
indices = resolve_genome_indices(
    genome_config='${genome_config}',
    genome_species='${genome_species}',
    genome_indices='${genome_indices_override}' if '${genome_indices_override}' != 'null' else None,
    download_indexes=${download_indexes}
)

# Write resolved indices to JSON
with open('resolved_indices.json', 'w') as f:
    json.dump(indices, f, indent=2)

print(f'Resolved indices for species: {list(indices.keys())}')
for species, index_path in indices.items():
    print(f'  {species}: {index_path}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > resolved_indices.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """
}
