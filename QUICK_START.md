# Quick start

This package contains the public MATLAB implementation and the saved outputs
for the revised MPHBS manuscript.

## Main algorithm

The final MPHBS implementation is:

```text
code/MPHBS_main.m
```

Function signature:

```matlab
[bestValue, bestPosition, convergenceCurve] = MPHBS_main( ...
    N, T, LB, UB, Dim, F_obj, K, Ni, rho_sub, w_best, gamma, ...
    fad_before_p2, fad_after_p2);
```

Final manuscript configurations:

```text
CEC2022 B-C1:
K = 7, Ni = 3, rho_sub = 0.25, w_best = 0.050,
gamma = 0.30, fad_before_p2 = true, fad_after_p2 = false

FIR A-C1:
K = 3, Ni = 9, rho_sub = 0.25, w_best = 0.050,
gamma = 0.30, fad_before_p2 = false, fad_after_p2 = true
```

## Related public implementations

```text
code/MPHBS_random_mirror.m   matched RANDOM mirror of MPHBS
code/MPHB_baseline.m         HBA-MPA backbone without the mediator
```

## Run the manuscript experiments

Start MATLAB in the `code/` directory:

```matlab
cd('/path/to/05_reproducibility_package_public_lean/code')
```

Run all six manuscript campaigns:

```matlab
run_all_manuscript_experiments
```

Run validation on the saved outputs:

```matlab
validate_campaign_outputs
validate_paired_seed_results
audit_final_fe_budgets
```

CSV-level validation can also be run from a terminal:

```bash
python3 code/validate_result_files.py
```

For complete details, see `README_REPRODUCIBILITY.md`.
