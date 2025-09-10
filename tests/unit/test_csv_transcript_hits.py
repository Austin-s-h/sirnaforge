import csv

import pytest

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.sirna import DesignParameters


@pytest.mark.unit
@pytest.mark.local_python
@pytest.mark.ci
def test_csv_contains_transcript_hit_columns(tmp_path):
    """Design from a FASTA with two identical transcripts, save CSV, and assert columns/values."""

    fasta = tmp_path / "test.fa"
    seq = "A" * 50
    fasta.write_text(f">tx1\n{seq}\n>tx2\n{seq}\n")

    params = DesignParameters(sirna_length=21, top_n=5)
    designer = SiRNADesigner(params)
    res = designer.design_from_file(str(fasta))

    out = tmp_path / "out.csv"
    res.save_csv(str(out))

    with out.open(newline="") as fh:
        reader = csv.DictReader(fh)
        # Ensure new columns are present
        assert "transcript_hit_count" in reader.fieldnames
        assert "transcript_hit_fraction" in reader.fieldnames

        rows = list(reader)
        # Rows should equal total_candidates
        assert len(rows) == res.total_candidates

        # At least one row should report hit_count == 2 (both transcripts identical)
        assert any(int(r.get("transcript_hit_count", 0)) == 2 for r in rows)
