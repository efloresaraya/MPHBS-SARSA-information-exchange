function run_all_manuscript_experiments
%RUN_ALL_MANUSCRIPT_EXPERIMENTS Execute the six paired-seed campaigns.

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
addpath(fullfile(codeDir, 'input_data'));
originalDir = pwd;
cleanup = onCleanup(@() cd(originalDir));
cd(codeDir);

run_CEC2022_sensitivity_common_seeds;
run_FIR_sensitivity_common_seeds;
run_CEC2022_ablation_BC1_common_seeds;
run_FIR_ablation_AC1_common_seeds;
run_CEC2022_external_BC1_common_seeds;
run_FIR_external_AC1_common_seeds;

export_official_block_statistics;
validate_paired_seed_results;
clear cleanup
end
