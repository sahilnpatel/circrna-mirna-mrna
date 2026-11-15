import sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from data_io import aligned_rows,read_matrix

def test_join_is_order_independent(tmp_path):
    from shutil import copyfile
    import csv
    for f in (ROOT/'data/source').glob('*.tsv'):copyfile(f,tmp_path/f.name)
    original=aligned_rows(tmp_path)
    path=tmp_path/'mcm5.tsv'
    rows=list(csv.reader(path.open(),delimiter='\t'))
    rows[0]=rows[0][::-1];rows[1]=[rows[1][0],*rows[1][1:][::-1]]
    with path.open('w') as s:csv.writer(s,delimiter='\t').writerows(rows)
    assert aligned_rows(tmp_path)==original


def test_duplicate_ids_fail(tmp_path):
    p=tmp_path/'bad.tsv';p.write_text('TCGA-AA-0001-01\tTCGA-AA-0001-01\nx\t1\t2\n')
    with pytest.raises(ValueError,match='Duplicate'):read_matrix(p)


def test_recovered_cohort():
    rows=aligned_rows(ROOT/'data/source')
    assert len(rows)==1129
    tumor=[r for r in rows if r['sample_type']=='01']
    assert len(tumor)==1017
    assert len({r['patient_id'] for r in tumor})==1017
    assert 'MCM5' in rows[0] and len([k for k in rows[0] if k.startswith(('GOBP_','REACTOME_'))])==11
    r=np.corrcoef([r['MCM5'] for r in rows],[r['mir885_activity_proxy'] for r in rows])[0,1]
    assert abs(r-(-.2809100973771341))<1e-10
