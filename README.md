# Avipoxvirus ELISA

Code and data for the FPV ELISA methods paper (Carlson & Forsman)

Semi-quantitative indirect ELISA for avipoxvirus antibodies in Florida Scrub-Jays, using a pooled logit-log standard curve to get ACU values from 4 plates of 30 birds each.

## What's here

**R scripts**
- `Supplementary_Code_S1.R` — the whole pipeline start to finish. Standard curve, mixed effects, ACU calculation, QC, exports. All OD values are hardcoded so you don't need the Excel files to run it.
- `fpv_figures.R` — main text + supplementary figures. Run after the main script.
- `fig_S1_4pl_comparison.R` — compares pooled logit-log to per-plate 4PL. Needs `minpack.lm`. Also run after the main script.

**Data**
- `acu_sample_means.csv` — the 120 birds, duplicate-averaged
- `acu_test_wells.csv` — all 240 individual wells
- `acu_qc_wells.csv` — QC triplicate data (high/mid/low × 4 plates)
- `standard_curve_data.csv` — 64 standard curve points

## Plate layout

Same for all 4 plates:

```
Col 1-2:  standard curve (128 -> 1 ACU, dup)
Col 3 AB: negative controls
Col 3 CE: QC high (FSJ-490), triplicate
Col 3 FH: QC mid (FSJ-475), triplicate
Col 4 AC: QC low (FSJ-491), triplicate
Col 4 D:  blank
Col 4 EH, Col 5-11: test samples, dup pairs A+B C+D E+F G+H
Col 12:   blanks
```

30 samples/plate × 4 plates = 120 birds.

## How to run

```r
source("Supplementary_Code_S1.R")   # everything
source("fpv_figures.R")              # figures (needs objects from above)
source("fig_S1_4pl_comparison.R")    # 4PL comparison (also needs above)
```

Packages install themselves if missing. Needs R ≥ 4.0.

## The short version of the method

Blank-subtract → normalize to relative binding (OD/OD_max) → logit transform → pool all 4 plates → OLS regression → invert to get ACU. Pooled equation: logit(y) = −3.345 + 2.461 × log10(ACU), R² = 0.957.

Mixed effects confirmed plate random effects are zero after normalization, so pooling is fine.

## Citation

Carlson, A.C. and Forsman, A.M. (2026). Development and validation of a semi-quantitative indirect ELISA for avipoxvirus antibody detection in a wild passerine using purified chorioallantoic membrane antigen.

