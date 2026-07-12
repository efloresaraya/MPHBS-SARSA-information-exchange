# MPHBS reproducibility package

This directory contains the MATLAB code and the complete saved outputs for
the six paired-seed campaigns reported in the revised MPHBS manuscript.

## Reproducibility status

- MATLAB release: R2024b
- Hardware used: Mac mini, Apple M4 Pro, 24 GB RAM
- Operating system: macOS 26.3.1
- Parallel pool: eight local process workers
- Runs: 30 paired runs per benchmark instance
- Population size: 30
- Master seed: `20260405`
- Proposed-method entry point: `code/MPHBS_main.m`
- CEC2022 objective calls: observed directly by the runner wrapper
- FIR objective calls: observed directly by the FIR objective counter

The same deterministic seed is used by every compared method within an
instance/run block:

```text
CEC2022 seed = 20260405 + 1000 * function_position + run
FIR seed     = 20260405 + 1000 * case_id + run
```

Seeds differ across runs and benchmark instances. This common-random-number
design supports run-level W/T/L comparisons and paired Wilcoxon signed-rank
tests within each function or case. Holm correction is applied across the
complete set of anchor-versus-competitor instance-level hypotheses in each
domain. Functions and cases, rather than repeated runs, are treated as the
independent blocks in the Friedman analysis. Kendall's W reports the omnibus
effect size, and anchor-based Friedman post hoc comparisons are Holm-adjusted
across competitors.

## Required MATLAB products

- MATLAB R2024b
- Parallel Computing Toolbox
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

The bundled official CEC2022 interface is
`code/cec22_test_func.mexmaca64`, which targets Apple Silicon macOS. The
official CEC2022 data files required for dimension 20 are provided under
`code/input_data/`. On another platform, compile or obtain the corresponding
official CEC2022 MEX binary and place it in `code/`.

## Experimental protocol

### CEC2022

- Functions: F1--F12
- Dimension: 20
- Bounds: `[-100, 100]`
- Function-evaluation budget: `10,000D = 200,000`
- Reported value: residual error from the official optimum
- Final MPHBS configuration: B-C1
- B-C1 parameters: `K=7`, `Ni=3`, FAD before the mediator,
  `w_best=0.05`, `gamma=0.30`, `rho_sub=0.25`

### FIR

- Cases: FIR1--FIR8
- Filter order: 30
- Dimension: 31
- Bounds: `[-1, 1]`
- FFT resolution: 2048
- Function-evaluation budget: 150,000
- Final MPHBS configuration: A-C1
- A-C1 parameters: `K=3`, `Ni=9`, FAD after the mediator,
  `w_best=0.05`, `gamma=0.30`, `rho_sub=0.25`

B-C1 and A-C1 are domain-specific configurations of the same MPHBS
architecture. They were identified in the internal paired-seed sensitivity
analysis and then fixed for the corresponding ablation and external
comparison.

The FE values stored in `task_results.csv` and the MAT workspaces are
observed objective-call counts, not nominal iteration-based estimates. The
budget is an upper limit. In CEC2022, all sensitivity, ablation, and external
methods consume between 199,614 and 200,000 FEs. In FIR, all methods consume
between 149,718 and 150,000 FEs. Thus, both domains use directly observed,
closely matched FE consumption without exceeding their declared budgets.

## Included campaigns

| Domain | Stage | Runner | Saved result directory |
|---|---|---|---|
| CEC2022 | Sensitivity | `run_CEC2022_sensitivity_common_seeds` | `CEC2022_sensitivity_common_seeds` |
| FIR | Sensitivity | `run_FIR_sensitivity_common_seeds` | `FIR_sensitivity_common_seeds` |
| CEC2022 | Ablation | `run_CEC2022_ablation_BC1_common_seeds` | `CEC2022_ablation_BC1_common_seeds` |
| FIR | Ablation | `run_FIR_ablation_AC1_common_seeds` | `FIR_ablation_AC1_common_seeds` |
| CEC2022 | External | `run_CEC2022_external_BC1_common_seeds` | `CEC2022_external_BC1_common_seeds` |
| FIR | External | `run_FIR_external_AC1_common_seeds` | `FIR_external_AC1_common_seeds` |

The six campaigns contain 19,800 completed optimization tasks. The external
comparisons include MPHB, HBA, MPA, EODE, RLTLBO, QQLMPA, GWO,
FDB-TLABC, and FTO.

The complete package was re-executed from a clean results directory on
2026-07-10. The regenerated outputs match the saved public outputs, excluding
wall-clock time by design. All six campaigns passed seed-pairing,
finite-output, FE-budget, anchor-identity, and block-level statistical
validation. See `documentation/VALIDATION_REPORT.md`.

## Directory structure

```text
05_reproducibility_package_public_lean/
|-- code/                 MATLAB algorithms, runners, validation scripts,
|                         CEC2022 MEX, and official input data
|-- results/              Complete outputs from all six campaigns
|-- manuscript_figures/   Final 400-dpi PNG and editable MATLAB FIG files
|-- documentation/        Package metadata and verification reports
|-- README_REPRODUCIBILITY.md
|-- FILE_MANIFEST_SHA256.txt
`-- PACKAGE_INVENTORY.txt
```

For a short guide to the main public algorithm entry point, see
`QUICK_START.md`.

Every result directory includes the available raw run-level CSV files,
summary CSV files, generated LaTeX tables, convergence figures, final
workspace, configuration snapshot, runner snapshot, and paired-statistical
exports. Intermediate checkpoint folders are not required to verify the
reported results and are stored separately in
`../05_optional_checkpoints_archive_public/`.

## Running the experiments

Start MATLAB in `code/`:

```matlab
cd('/path/to/05_reproducibility_package_public_lean/code')
```

Run one campaign by invoking its runner without the `.m` extension. Run all
six campaigns with:

```matlab
run_all_manuscript_experiments
```

The runners write to `../results/` and create checkpoints during execution.
After the six campaigns, `run_all_manuscript_experiments` automatically
invokes `export_official_block_statistics`. To perform a completely fresh
execution, first copy the package and remove or rename the corresponding
result directory in that copy.

## Validating the saved outputs

The complete MATLAB validation checks successful finite outcomes, dimensions,
observed FE budgets, the common seed schedule, and exact anchor equality
between ablation and external experiments:

```matlab
validate_paired_seed_results
```

Static code and dependency checks:

```matlab
validate_campaign_outputs
```

CSV-level validation can also be run from a terminal:

```bash
python3 code/validate_result_files.py
```

Recompute all official statistical outputs directly from the raw run-level
CSV files:

```matlab
export_official_block_statistics
```

Audit every observed FE count against the declared CEC2022 and FIR budgets,
the full-precision MAT arrays, and a minimum 99.5% utilization threshold:

```matlab
audit_final_fe_budgets
```

The saved report is
`documentation/FE_BUDGET_AUDIT.txt`.

The official exporter verifies 30 unique common-seed pairs for every
anchor-competitor function or case before computing a signed-rank test.

The public entry point is `MPHBS_main.m`. MATLAB parse/dependency checks and
paired-result validations are provided by `validate_campaign_outputs.m` and
`validate_paired_seed_results.m`.

## Statistical outputs

Each campaign contains `statistics_official_block_level/` with:

- paired Wilcoxon signed-rank tests within every function or case;
- domain-wide Holm-adjusted instance-level p-values;
- function/case mean-based W/T/L summaries;
- Friedman average ranks and anchor-based post hoc comparisons;
- Holm-adjusted post hoc p-values;
- Friedman omnibus p-values and Kendall's W; and
- a MAT snapshot of the complete statistical export.

The official exporter removes runner-native pooled inferential summaries
because repeated runs are not independent benchmark problems. Historical
copies are retained outside the public package in the project provenance
archive.

No independent-sample test is used by the public paired-seed runners.

## Integrity

`FILE_MANIFEST_SHA256.txt` records the SHA-256 checksum of every package file
except the manifest itself. `PACKAGE_INVENTORY.txt` summarizes file counts,
sizes, critical artifacts, and validation results.
