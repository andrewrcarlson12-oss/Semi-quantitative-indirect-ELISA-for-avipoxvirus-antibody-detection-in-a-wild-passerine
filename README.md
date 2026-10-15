# Avipoxvirus ELISA

Code and data for the FPV ELISA methods paper (Carlson & Forsman).

Semi-quantitative indirect ELISA for avipoxvirus antibodies in Florida Scrub-Jays. Pooled logit-log standard curve, ACU values back-calculated across 4 plates of 30 birds each.

## What's here

**Code**
- `Supplementary Code S1`: the whole pipeline. Standard curve fit, mixed-effects sanity checks, ACU back-calculation, QC stats, exports. All OD values are hardcoded from the Gen5 plate reader exports, so you don't need any Excel files to run it.
- `Supplementary Code S2`: the experiments that ran on other plates. Antigen depletion test, dilutional parallelism, pre-absorption, spike and recovery, urea avidity, lot-to-lot comparison. ODs are hardcoded the same way. It prints the results for each one and doesn't write any files.

**Data**
- `ACU Sample Values`: 120 birds, duplicate-averaged ACU and mean OD
- `ACU Test Wells`: all 240 individual wells with raw OD, blank-subtracted OD, and back-calculated ACU
- `ACU QC Wells`: QC triplicate data (high/mid/low × 4 plates)
- `Standard Curve`: 64 standard curve points (8 dilutions × 2 reps × 4 plates), raw OD through logit-transformed relative binding

Running S1 regenerates the four data files. The values are identical to what's already in the repo.

## Plate layout

Same on all 4 plates:

```
Col 1-2:  standard curve (128 → 1 ACU, duplicate)
Col 3 AB: negative controls
Col 3 CE: QC high (FSJ-490), triplicate
Col 3 FH: QC mid (FSJ-475), triplicate
Col 4 AC: QC low (FSJ-491), triplicate
Col 4 D:  blank
Col 4 EH, Col 5-11: test samples, duplicate pairs A+B, C+D, E+F, G+H
Col 12:   blanks
```

30 test samples per plate × 4 plates = 120 birds.

## How to run

```r
source("Supplementary Code S1")
source("Supplementary Code S2")
```

That's it. S1 needs `lme4` and `dplyr`. It runs the OLS fit, runs the mixed-effects model to confirm plate random effects are negligible, computes intra- and inter-assay CVs, flags the Grubbs outlier, and writes the four CSVs. S2 only needs base R, but its lot comparison reads `Standard Curve`, so run it from the repo folder. Needs R ≥ 4.0.

## Method in one paragraph

Blank-subtract each well, normalize to relative binding (`bs_OD / OD_max`), logit-transform, pool all four plates, fit OLS, invert to get ACU.

Pooled equation: `logit(y) = −3.345 + 2.461 × log10(ACU)`, R² = 0.957.

The mixed-effects model exists to confirm what you can already see from the per-plate slopes (CV = 4.1%): plate random effects collapse to zero after normalization, so pooling is fine. The boundary singular fit warning is the point, not a problem.

## Citation

Carlson AR, Forsman AM. 2026. Development and analytical validation of a semi-quantitative indirect ELISA for detection of anti-avipoxvirus antibodies in Florida Scrub-Jays (*Aphelocoma coerulescens*) using purified chorioallantoic membrane antigen. *Conservation Physiology*.
