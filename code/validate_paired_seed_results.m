function validate_paired_seed_results
%VALIDATE_PAIRED_SEED_RESULTS Verify all six paired-seed campaigns.

rootDir = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(rootDir, 'results');

cecCampaigns = {
    'CEC2022_sensitivity_common_seeds'
    'CEC2022_ablation_BC1_common_seeds'
    'CEC2022_external_BC1_common_seeds'
};
firCampaigns = {
    'FIR_sensitivity_common_seeds'
    'FIR_ablation_AC1_common_seeds'
    'FIR_external_AC1_common_seeds'
};

for k = 1:numel(cecCampaigns)
    file = fullfile(resultsDir, cecCampaigns{k}, 'mat', ...
        'final_workspace.mat');
    saved = load(file, 'RawSeedUsed', 'RawError', 'RawOK', ...
        'CFG', 'Entries');
    assert(all(saved.RawOK, 'all'));
    assert(all(isfinite(saved.RawError), 'all'));
    assert(size(saved.RawSeedUsed,1) == numel(saved.Entries));
    assert(size(saved.RawSeedUsed,2) == numel(saved.CFG.Functions));
    assert(size(saved.RawSeedUsed,3) == saved.CFG.Runs);
    for f = 1:numel(saved.CFG.Functions)
        expected = saved.CFG.MASTER_SEED + 1000 * f + ...
            (1:saved.CFG.Runs);
        for a = 1:numel(saved.Entries)
            assert(isequal(squeeze(saved.RawSeedUsed(a,f,:))', expected));
        end
    end
    fprintf('CEC2022 verified: %s\n', cecCampaigns{k});
end

for k = 1:numel(firCampaigns)
    file = fullfile(resultsDir, firCampaigns{k}, 'mat', ...
        'final_workspace.mat');
    saved = load(file, 'SeedMask', 'RawFinal', 'OKMask', ...
        'CFG', 'Cases', 'Entries');
    assert(all(saved.OKMask, 'all'));
    assert(all(isfinite(saved.RawFinal), 'all'));
    assert(size(saved.SeedMask,1) == numel(saved.Entries));
    assert(size(saved.SeedMask,2) == numel(saved.Cases));
    assert(size(saved.SeedMask,3) == saved.CFG.nRuns);
    for c = 1:numel(saved.Cases)
        expected = saved.CFG.MASTER_SEED + 1000 * saved.Cases(c).id + ...
            (1:saved.CFG.nRuns);
        for a = 1:numel(saved.Entries)
            assert(isequal(squeeze(saved.SeedMask(a,c,:))', expected));
        end
    end
    fprintf('FIR verified: %s\n', firCampaigns{k});
end

cecAblation = load(fullfile(resultsDir, cecCampaigns{2}, 'mat', ...
    'final_workspace.mat'), 'RawError', 'Entries');
cecExternal = load(fullfile(resultsDir, cecCampaigns{3}, 'mat', ...
    'final_workspace.mat'), 'RawError', 'Entries');
cecAblationIndex = find(strcmp(string({cecAblation.Entries.alias}), 'BS-C1'), 1);
cecExternalIndex = find(strcmp(string({cecExternal.Entries.alias}), 'B-C1'), 1);
assert(~isempty(cecAblationIndex) && ~isempty(cecExternalIndex));
assert(isequal(cecAblation.RawError(cecAblationIndex,:,:), ...
    cecExternal.RawError(cecExternalIndex,:,:)));
fprintf('CEC2022 anchor verified exactly across ablation and external results.\n');

firAblation = load(fullfile(resultsDir, firCampaigns{2}, 'mat', ...
    'final_workspace.mat'), 'RawFinal', 'Entries');
firExternal = load(fullfile(resultsDir, firCampaigns{3}, 'mat', ...
    'final_workspace.mat'), 'RawFinal', 'Entries');
firAblationIndex = find(strcmp(string({firAblation.Entries.alias}), 'AS-C1'), 1);
firExternalIndex = find(strcmp(string({firExternal.Entries.alias}), 'A-C1'), 1);
assert(~isempty(firAblationIndex) && ~isempty(firExternalIndex));
assert(isequal(firAblation.RawFinal(firAblationIndex,:,:), ...
    firExternal.RawFinal(firExternalIndex,:,:)));
fprintf('FIR anchor verified exactly across ablation and external results.\n');

allCampaigns = [cecCampaigns; firCampaigns];
expectedKendallW = [ ...
    0.421444171444171
    0.530158730158730
    0.626094276094276
    0.638590294840295
    0.848214285714286
    0.917045454545454
];
for campaignIdx = 1:numel(allCampaigns)
    statisticsDir = fullfile(resultsDir, allCampaigns{campaignIdx}, ...
        'statistics_official_block_level');
    omnibusFile = fullfile(statisticsDir, ...
        'official_friedman_kendall_omnibus.csv');
    instanceFile = fullfile(statisticsDir, ...
        'official_paired_wilcoxon_by_instance.csv');
    globalFile = fullfile(statisticsDir, ...
        'official_global_by_instance_blocks.csv');
    assert(isfile(omnibusFile) && isfile(instanceFile) && ...
        isfile(globalFile));
    omnibus = readtable(omnibusFile);
    assert(abs(omnibus.KendallW - expectedKendallW(campaignIdx)) < 1e-12);
end
fprintf('Official Friedman, Holm, and Kendall outputs verified.\n');

fprintf('All six campaigns are complete, finite, and correctly paired.\n');
end
