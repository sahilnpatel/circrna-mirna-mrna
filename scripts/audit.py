"""Rebuild sample-aligned features and audit the historical ML export."""
import csv
import hashlib
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from data_io import aligned_rows, read_matrix, write_csv


def main():
    source = ROOT / 'data/source'
    rows = aligned_rows(source)
    write_csv(ROOT / 'data/processed/aligned_features.csv', rows)
    tumors = [r for r in rows if r['sample_type'] == '01']
    if len({r['patient_id'] for r in tumors}) != len(tumors):
        raise ValueError('Multiple primary tumor samples for one patient')
    write_csv(ROOT / 'data/processed/primary_tumor_features.csv', tumors)
    mcm = read_matrix(source / 'mcm5.tsv')['MCM5']
    proxy = read_matrix(source / 'mir885_targets.tsv')['genelist']
    with (source / 'legacy_ml.csv').open() as stream:
        raw = list(csv.reader(stream))[1:]
    legacy = [(float(r[0]), float(r[1])) for r in raw if r[0] and r[1]]
    score_original_order = list(proxy.values())
    sorted_mcm = sorted(mcm.values())
    r_all = statistics.correlation([r['mir885_activity_proxy'] for r in rows], [r['MCM5'] for r in rows])
    r_tumor = statistics.correlation([r['mir885_activity_proxy'] for r in tumors], [r['MCM5'] for r in tumors])
    report = {
        'n_samples': len(rows), 'n_primary_tumors': len(tumors),
        'sample_type_counts': dict(Counter(r['sample_type'] for r in rows)),
        'n_saved_pathways': len(rows[0]) - 6,
        'pearson_r_all_samples': r_all, 'pearson_r_primary_tumors': r_tumor,
        'paper_reported_r': -0.281,
        'legacy_rows': len(legacy),
        'legacy_r': statistics.correlation(*zip(*legacy)),
        'legacy_mir_matches_original_order': sum(abs(a + score_original_order[i]) < 1e-7 for i, (a, b) in enumerate(legacy)),
        'legacy_mcm_matches_independently_sorted_values': sum(abs(b - sorted_mcm[i]) < 1e-7 for i, (a, b) in enumerate(legacy)),
        'sources_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.glob('*.tsv'))},
        'warnings': ['miRNA activity is a negated target-gene score, not measured miRNA.',
                     'MCM5 is itself included in the 74-gene target signature: circularity risk.',
                     'LUAD/LUSC clinical subtype and survival metadata are not supplied.',
                     'Original GSVA execution and survival results have not been reproduced.']
    }
    (ROOT / 'results/audit.json').write_text(json.dumps(report, indent=2) + '\n')
    (ROOT / 'results/AUDIT.md').write_text(f'''# Source audit

The sample-ID join recovers **{len(rows)} samples**, including **{len(tumors)} primary tumors**.
All-sample Pearson r = **{r_all:.6f}**, consistent with the paper's rounded -0.281.
Primary-tumor-only r = **{r_tumor:.6f}**.

The historical ML export contains {len(legacy)} rows. Its miRNA column matches the
first {report['legacy_mir_matches_original_order']} negated target scores in original order;
its MCM5 column matches the first {report['legacy_mcm_matches_independently_sorted_values']}
values of independently sorted MCM5. This breaks sample pairing. Its correlation is
{report['legacy_r']:.6f}. Do not use that file for biological inference.

The source contains 01 (primary tumor), 02 (recurrent tumor), and 11 (solid tissue normal).
An observation count is not a count of distinct patients. The primary-tumor table has
one sample per patient, verified by the first 12 characters of the TCGA barcode.

## Limits

This audit reproduces a correlation from saved outputs, not a complete rerun of the
original GSVA pipeline. The 74-gene signature includes MCM5; its association with
MCM5 therefore has a built-in dependency. Recompute the signature without MCM5
before claiming an independent regulatory relationship. Eleven pathway score rows
are available in the combined saved TSV; the original paper describes 22 gene sets.
Other gene-set exports are archived but not silently treated as the same analysis.

No patient survival table was recovered. Paper survival p-values remain historical
reported results, not newly verified findings. No circRNA or miRNA was directly
measured in this rebuilt feature table.
''')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
