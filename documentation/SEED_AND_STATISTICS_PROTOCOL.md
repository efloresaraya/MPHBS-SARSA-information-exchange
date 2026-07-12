# Seed and statistical protocol

All methods use the same seed within each benchmark-instance/run block.
The random-number generator is MATLAB `twister`.

```text
CEC2022: seed = 20260405 + 1000 * function_position + run
FIR:     seed = 20260405 + 1000 * case_id + run
```

The run index ranges from 1 to 30. Algorithms are paired by benchmark
instance, run index, and seed. Paired runs are compared within an instance,
whereas functions or cases are the independent blocks for domain-level
inference. The analysis reports:

- arithmetic mean and standard deviation by benchmark instance;
- average ranks computed from instance-level means;
- run-level wins, ties, and losses against the declared anchor;
- paired Wilcoxon signed-rank tests within each function or case;
- domain-wide Holm correction across all anchor-competitor instance tests;
- Friedman omnibus tests over function/case mean ranks;
- Kendall's W as the omnibus effect size;
- anchor-based Friedman post hoc comparisons with Holm correction across
  competitors; and
- convergence curves averaged over the 30 paired runs.

CEC2022 uses B-C1 as the final anchor. FIR uses A-C1 as the final anchor.
The corresponding ablation aliases are BS-C1 and AS-C1.

The pairing unit is the seed-labelled replicate for a fixed benchmark
instance. Different algorithms can consume different numbers of random
variates after initialization; pairing therefore denotes a matched
replicate, not an assertion that every internal random draw is identical.
Repeated runs are not pooled across heterogeneous functions or cases for
domain-level inference.

The complete analysis is reproduced by:

```matlab
export_official_block_statistics
```

Its outputs are stored under
`results/<campaign>/statistics_official_block_level/`.
