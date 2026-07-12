campaignCode = fileparts(mfilename('fullpath'));
addpath(campaignCode);

runnerFiles = {
    'run_CEC2022_sensitivity_common_seeds.m'
    'run_FIR_sensitivity_common_seeds.m'
    'run_CEC2022_ablation_BC1_common_seeds.m'
    'run_FIR_ablation_AC1_common_seeds.m'
    'run_CEC2022_external_BC1_common_seeds.m'
    'run_FIR_external_AC1_common_seeds.m'
    'export_official_block_statistics.m'
    'run_all_manuscript_experiments.m'
    'validate_paired_seed_results.m'
    'audit_final_fe_budgets.m'
    };

for i = 1:numel(runnerFiles)
    filePath = fullfile(campaignCode, runnerFiles{i});
    issues = checkcode(filePath, '-id');
    fprintf('%s: parsed successfully (%d Code Analyzer notices).\n', ...
        runnerFiles{i}, numel(issues));
end

assert(strcmp(which('MPHBS_main'), fullfile(campaignCode, 'MPHBS_main.m')));
assert(~isempty(which('cec22_test_func')));
fprintf('Campaign dependencies resolve from: %s\n', campaignCode);
