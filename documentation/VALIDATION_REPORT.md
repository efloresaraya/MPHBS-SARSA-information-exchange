# Validation report - 2026-07-10

## Official campaigns

The lean public reproducibility package contains six validated common-seed campaigns:

- CEC2022 sensitivity, B-C1 ablation, and B-C1 external comparison;
- FIR sensitivity, A-C1 ablation, and A-C1 external comparison.

Together they contain 19,800 successful optimization tasks. Every
benchmark-instance/run block uses one common seed across all compared
methods. All saved objective values are finite. Intermediate checkpoint
folders are archived separately because they are not required to verify the
reported CSV files, final workspaces, figures, or statistical exports.

## Function evaluations

CEC2022 uses a maximum budget of 200,000 directly observed objective calls.
The final saved campaigns contain 199,614--200,000 calls per task. FIR uses
a maximum budget of 150,000 calls and contains 149,718--150,000 calls per
task. No validated task exceeds its declared budget.

## Final configurations

- CEC2022: B-C1 (`K=7`, `Ni=3`, BEFORE, `w_best=0.05`,
  `gamma=0.30`, `rho_sub=0.25`);
- FIR: A-C1 (`K=3`, `Ni=9`, AFTER, `w_best=0.05`,
  `gamma=0.30`, `rho_sub=0.25`).

B-C1 is exactly identical across the corresponding CEC2022 ablation and
external campaigns. A-C1 is exactly identical across the FIR ablation and
external campaigns.

## Statistical validation

`export_official_block_statistics.m` reconstructs the official analysis from
the raw `task_results.csv` files. It verifies 30 unique common-seed pairs for
every anchor-competitor function or case and exports:

- paired Wilcoxon tests within each instance;
- domain-wide Holm correction;
- function/case mean-based W/T/L and average ranks;
- Friedman omnibus and anchor-based post hoc tests;
- Holm correction across post hoc competitors; and
- Kendall's W.

The reproduced omnibus values are:

| Campaign | Friedman p | Kendall's W |
|---|---:|---:|
| CEC2022 sensitivity | 1.13672e-10 | 0.421444 |
| CEC2022 ablation | 6.48017e-06 | 0.530159 |
| CEC2022 external | 4.45440e-11 | 0.626094 |
| FIR sensitivity | 8.07202e-11 | 0.638590 |
| FIR ablation | 2.46025e-06 | 0.848214 |
| FIR external | 9.10310e-11 | 0.917045 |

These values reproduce the revised manuscript. Older pooled summaries are
retained only in the project provenance archive outside the public package
and are not used for inference.

## Figures

The four manuscript convergence figures correspond to the official B-C1 and
A-C1 campaigns. PNG files are exported at 400 dpi, and editable MATLAB FIG
files are included.

## Independent Python cross-check

A separate Python audit reconstructed the paired seeds, tolerant W/T/L counts, approximate Wilcoxon signed-rank p-values, and Holm-adjusted p-values from the raw `task_results.csv` files. The maximum absolute difference from the official CSV exports was `5.829e-16` for raw p-values and `1.887e-15` for Holm-adjusted p-values across 600 hypotheses.
