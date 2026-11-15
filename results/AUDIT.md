# Source audit

The sample-ID join recovers **1129 samples**, including **1017 primary tumors**.
All-sample Pearson r = **-0.280910**, consistent with the paper's rounded -0.281.
Primary-tumor-only r = **-0.248850**.

The historical ML export contains 1000 rows. Its miRNA column matches the
first 1000 negated target scores in original order;
its MCM5 column matches the first 1000
values of independently sorted MCM5. This breaks sample pairing. Its correlation is
0.027638. Do not use that file for biological inference.

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
