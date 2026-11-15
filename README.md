# circRNA–miRNA–mRNA regulation in non-small cell lung cancer

This repository accompanies my paper, **[Unraveling a CircRNA–miRNA–mRNA Axis:
A Potential Therapeutic Target for Non-small Cell Lung Cancer](archive/paper/Research_Paper.pdf)**.
I began the study the summer after my freshman year of high school. It examines
the proposed **circFNDC3B–miR-885-3p–MCM5** axis and its relationship to NSCLC.

[Paper](archive/paper/Research_Paper.pdf) ·
[Research presentation](archive/paper/Research_Presentation.pdf) ·
[Manuscript listing student and mentor](archive/paper/Patel_Subramanian.pdf)

## Research question

Could a circRNA–miRNA–mRNA regulatory axis help identify a potential therapeutic
target in NSCLC? The paper follows circFNDC3B, miR-885-3p, and the target gene MCM5,
then examines pathways associated with MCM5 expression.

## Approach and findings in the paper

CircNet was used to identify candidate interactions. GEPIA2 analyses examined
survival associations in LUAD and LUSC, while TCGA lung expression data and GSVA
were used to investigate target-gene and pathway enrichment. Samples were also
compared by high and low MCM5 expression.

The paper reports a negative association between the inferred miR-885-3p activity
and MCM5 (r = −0.281), alongside survival associations and differences involving
Ras protein signal transduction and IGF transport and uptake by IGF-binding
proteins. These computational findings support a candidate axis for further study;
they do not establish a causal mechanism or a validated treatment target.

## What is included here

The original paper, presentation, R scripts, gene sets, and saved outputs are in
`archive/` and `data/source/`. A newer Python analysis checks sample matching,
reproduces the saved-output correlation, and reexamines the available pathway
scores. The original study is complete; reproducing every step from raw data is
still in progress. The [reproduction notes](docs/REPRODUCTION_STATUS.md) distinguish
results reported in the paper from results rerun in this repository.

## Start here

1. Read [the source audit](results/AUDIT.md).
2. Review [what is reproduced and what is not](docs/REPRODUCTION_STATUS.md).
3. Run the commands below, then compare the generated figures with the paper.
4. Continue into the linked [ML extension](https://github.com/sahilnpatel/nsclc-rna-ml).

## What the saved data show

- 1,129 samples: 1,017 primary tumor, 110 solid-tissue normal, two recurrent tumor.
- Correctly joined all-sample Pearson correlation: **−0.280910**, matching the paper's
  rounded **−0.281**. Primary-tumor-only correlation: **−0.248850**.
- The later ML spreadsheet independently sorted MCM5 and lost sample matching.
  Its 1,000 rows are not a valid biological feature table.
- The inferred miRNA activity is a negated target-gene enrichment score. It is not
  measured miRNA expression. Its 74-gene signature includes MCM5, creating a
  dependency that must be addressed in follow-up work.

![Sample-aligned study data](results/all_samples_correlation.svg)

Blue: primary tumor. Orange: normal tissue. Green: recurrent tumor.

## Run

Use Python 3.12. From a terminal:

```bash
git clone https://github.com/sahilnpatel/circrna-mirna-mrna.git
cd circrna-mirna-mrna
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate with `.venv\Scripts\activate` instead. See
[environment notes](docs/ENVIRONMENT.md) for the tested versions and macOS settings.

The source audit needs only Python's standard library:

```bash
python3 scripts/audit.py
```

For statistical comparisons and figures, follow [environment setup](docs/ENVIRONMENT.md), then:

```bash
python scripts/analyze_saved_outputs.py
python -m pytest tests -q
```

Optional modern GSVA rerun, after installing R and Bioconductor GSVA:

```bash
Rscript scripts/recompute_gsva.R /path/to/TCGA.LUNG.sampleMap_HiSeqV2.txt recomputed
```

The R script is prepared but has not been run. It uses `gsvaParam` and records
gene-set coverage and package versions. It also recomputes the target signature
without MCM5, an important follow-up check. The original GSVA package versions
and browser-based GEPIA2 survival analyses have not been reproduced.

## Reproducing the saved analyses

`audit.py` joins the original exports by sample ID and writes the aligned tables
and an audit report. `analyze_saved_outputs.py` calculates correlations for the
full cohort and primary tumors separately. It then divides primary tumors at the
median MCM5 value and compares 11 saved pathway scores using Welch tests with
Benjamini–Hochberg correction.

These pathway comparisons are a new descriptive reanalysis. They are not a full
reproduction of the paper's 22-gene-set analysis or its GEPIA2 survival plots.
Tables, figures, and machine-readable summaries are written to `results/`.
The [validation record](docs/VALIDATION.md) links to the successful GitHub run.

## Files

| Folder | Contents |
|---|---|
| `archive/code` | Original R scripts, including historical errors |
| `archive/paper` | Supplied paper and presentation versions |
| `archive/gene_sets` | Original gene-set files and additional score workbooks |
| `data/source` | Small source exports required to reproduce the audit |
| `data/processed` | Sample-ID-aligned tables generated by the audit |
| `scripts` / `src` | New reproducible analysis and strict readers |
| `results` | Regenerated tables, audit, and figures |
| `docs` | Provenance, scope, limitations, and learning guidance |
| `tests` | Sample alignment and statistical regression checks |

See [source hashes](docs/source_manifest.json), [credit and rights](docs/RIGHTS_AND_CREDIT.md),
and [references](docs/REFERENCES.md). Correlation and computational prediction do
not establish a causal mechanism or clinical utility.
