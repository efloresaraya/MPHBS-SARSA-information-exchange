from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CAMPAIGNS = {
    "cec_sensitivity": (
        "CEC2022_sensitivity_common_seeds",
        6120,
        17,
        "FunctionID",
        "Error",
        {"A-C1", "A-C2", "A-C3", "A-C4",
         "B-C1", "B-C2", "B-C3", "B-C4", "MPHB",
         "RA-C1", "RA-C2", "RA-C3", "RA-C4",
         "RB-C1", "RB-C2", "RB-C3", "RB-C4"},
        0.421444171444171,
    ),
    "fir_sensitivity": (
        "FIR_sensitivity_common_seeds",
        4080,
        17,
        "CaseID",
        "BestFit",
        {"A-C1", "A-C2", "A-C3", "A-C4",
         "B-C1", "B-C2", "B-C3", "B-C4", "MPHB",
         "RA-C1", "RA-C2", "RA-C3", "RA-C4",
         "RB-C1", "RB-C2", "RB-C3", "RB-C4"},
        0.638590294840295,
    ),
    "cec_ablation": (
        "CEC2022_ablation_BC1_common_seeds",
        2160,
        6,
        "FunctionID",
        "Error",
        {"BS-C1", "BR-C1", "R-best", "MPHB", "HBA", "MPA"},
        0.530158730158730,
    ),
    "fir_ablation": (
        "FIR_ablation_AC1_common_seeds",
        1440,
        6,
        "CaseID",
        "BestFit",
        {"AS-C1", "AR-C1", "R-best", "MPHB", "HBA", "MPA"},
        0.848214285714286,
    ),
    "cec_external": (
        "CEC2022_external_BC1_common_seeds",
        3600,
        10,
        "FunctionID",
        "Error",
        {"B-C1", "MPHB", "HBA", "MPA", "EODE", "RLTLBO",
         "QQLMPA", "GWO", "FDB-TLABC", "FTO"},
        0.626094276094276,
    ),
    "fir_external": (
        "FIR_external_AC1_common_seeds",
        2400,
        10,
        "CaseID",
        "BestFit",
        {"A-C1", "MPHB", "HBA", "MPA", "EODE", "RLTLBO",
         "QQLMPA", "GWO", "FDB-TLABC", "FTO"},
        0.917045454545454,
    ),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


loaded: dict[str, list[dict[str, str]]] = {}
for key, (
    dirname,
    expected_rows,
    method_count,
    instance_col,
    value_col,
    aliases,
    expected_kendall_w,
) in CAMPAIGNS.items():
    campaign = ROOT / "results" / dirname
    task_path = campaign / "csv" / "task_results.csv"
    summary_path = campaign / "csv" / "summary_global.csv"
    workspace_path = campaign / "mat" / "final_workspace.mat"
    statistics_dir = campaign / "statistics_official_block_level"

    assert task_path.is_file(), task_path
    assert summary_path.is_file(), summary_path
    assert workspace_path.is_file(), workspace_path
    instance_path = statistics_dir / "official_paired_wilcoxon_by_instance.csv"
    global_path = statistics_dir / "official_global_by_instance_blocks.csv"
    omnibus_path = statistics_dir / "official_friedman_kendall_omnibus.csv"
    assert instance_path.is_file(), instance_path
    assert global_path.is_file(), global_path
    assert omnibus_path.is_file(), omnibus_path

    omnibus = read_csv(omnibus_path)
    assert len(omnibus) == 1
    kendall_w = float(omnibus[0]["KendallW"])
    assert math.isclose(kendall_w, expected_kendall_w, abs_tol=1e-12)

    rows = read_csv(task_path)
    loaded[key] = rows
    assert len(rows) == expected_rows, (key, len(rows), expected_rows)
    assert {row["Alias"] for row in rows} == aliases
    assert all(row["OK"] in {"1", "true", "True"} for row in rows)
    assert all(math.isfinite(float(row[value_col])) for row in rows)

    blocks: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        blocks[(row[instance_col], row["RunIdx"])].append(row)
    assert all(len(block) == method_count for block in blocks.values())
    assert all(len({row["Seed"] for row in block}) == 1 for block in blocks.values())

    summary = read_csv(summary_path)
    assert len(summary) == method_count
    print(f"{key}: PASS ({len(rows)} tasks, {len(blocks)} paired blocks)")


def anchor_map(rows: list[dict[str, str]], alias: str, instance_col: str, value_col: str):
    return {
        (row[instance_col], row["RunIdx"], row["Seed"]): float(row[value_col])
        for row in rows
        if row["Alias"] == alias
    }


cec_ablation = anchor_map(loaded["cec_ablation"], "BS-C1", "FunctionID", "Error")
cec_external = anchor_map(loaded["cec_external"], "B-C1", "FunctionID", "Error")
assert cec_ablation == cec_external
print("CEC anchor reproducibility: PASS (B-C1 exact match across 360 paired runs)")

fir_ablation = anchor_map(loaded["fir_ablation"], "AS-C1", "CaseID", "BestFit")
fir_external = anchor_map(loaded["fir_external"], "A-C1", "CaseID", "BestFit")
assert fir_ablation == fir_external
print("FIR anchor reproducibility: PASS (A-C1 exact match across 240 paired runs)")
