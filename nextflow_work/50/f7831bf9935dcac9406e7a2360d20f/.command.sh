#!/bin/bash -euo pipefail
python3 -c "
import sys
sys.path.insert(0, '/home/runner/work/sirnaforge/sirnaforge/src/sirnaforge/pipeline/nextflow/workflows/../src')
from sirnaforge.core.off_target import validate_and_write_sequences

# Validate siRNA candidates
total, valid, errors = validate_and_write_sequences(
    input_file='input_candidates.fasta',
    output_file='validated_candidates.fasta',
    expected_length=21
)

# Write validation report
with open('validation_report.txt', 'w') as f:
    f.write(f'Total candidates: {total}\n')
    f.write(f'Valid candidates: {valid}\n')
    f.write(f'Invalid candidates: {total - valid}\n')
    if errors:
        f.write('\nErrors:\n')
        for error in errors:
            f.write(f'  {error}\n')

print(f'Validated {valid} out of {total} candidates')
"

    cat <<-END_VERSIONS > versions.yml
    "SIRNAFORGE_OFFTARGET:SIRNA_OFFTARGET_ANALYSIS:PREPARE_CANDIDATES":
        python: $(python --version | sed 's/Python //g')
        biopython: $(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
