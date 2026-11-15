# Reproduction status

| Component | Status | Evidence / next requirement |
|---|---|---|
| Original study papers and presentation | Archived | `archive/paper` |
| Original R code | Archived | Three scripts; preserved with historical errors |
| miRNA target signature | Recovered | 74 genes extracted from original R script |
| TCGA expression matrix | Located, not committed | 167 MB; SHA-256 in manifest |
| R workspace | Located, not committed | 65 MB; SHA-256 in manifest |
| Sample alignment | Rebuilt and tested | Join on normalized TCGA sample IDs |
| Reported correlation −0.281 | Reproduced from saved outputs | `results/audit.json` |
| Normal/tumor composition | Audited | Barcode sample codes: 01 / 02 / 11 |
| Pathway comparisons | New descriptive reanalysis | 11 aligned saved pathways, Welch tests + BH adjustment |
| Full original 22-gene-set analysis | Incomplete reproduction | Additional workbooks/GMTs archived; map each original set and version |
| Original GSVA computation | Not rerun | R unavailable; original package version absent |
| Modern GSVA and leave-MCM5-out signature | Script prepared, not executed | `scripts/recompute_gsva.R` |
| CircNet ranking and interaction confidence | Historical paper evidence only | Need original database version/export |
| GEPIA2 survival plots/p-values | Historical paper evidence only | Need clinical outcomes and exact cutoffs/cohorts |
| Independent biological validation | Not performed | Direct assays or independent cohort needed |

## Historical code audit

`Data_Set_Code.R` reads a relative expression matrix. `GSVAmRNAfunction.R` uses an
interactive file chooser and the old GSVA interface. The miRNA script contains
`rite.table` (a typo) and a VEGF export unrelated to the target analysis. These
files are source history, not the supported execution path.

The saved miRNA TSV contains one unlabeled duplicate numeric row after blank
rows. The reader discards it only after verifying exact equality to an already
read feature. Other malformed rows fail rather than being silently ignored.
The combined pathway table also contains MCM5: it is verified against the
original MCM5 export and excluded from the pathway feature list.

## First hands-on task

Run `python scripts/audit.py`. Explain why matching sample IDs is different from
matching row numbers. Compare the first row of the old ML sheet with the sample
ID and MCM5 value in the aligned table. Then explain why 1,129 samples cannot
simply be described as 1,129 independent tumor patients.
