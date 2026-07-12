function run_CEC2022_ablation_BC1_common_seeds
% =========================================================================
% RUN_CEC2022_SENSITIVITY_COMMON_SEEDS
% =========================================================================
% CEC2022 ablation for the best-ranked SARSA sensitivity setting (B-C1).
%
% Methods
% -------
%   SARSA  : B-C1..B-C4, A-C1..A-C4
%   RANDOM : RB-C1..RB-C4, RA-C1..RA-C4
%   BASE   : MPHB
%
% Outputs
% -------
%   CSV/TEX:
%     cec_alias_map
%     cec_sarsa_sensitivity
%     cec_random_sensitivity
%     cec_anchor_stats
%
%   Figures:
%     cec_convergence_panels_logFE.png / .fig
%
% Design principles
% -----------------
% - Rank for each F is computed from the mean errors of the FULL method set.
% - AvgRank is the mean of the per-function global ranks.
% - Panel tables (SARSA / RANDOM) only FILTER displayed columns; they do NOT
%   recompute ranks on the displayed subset.
% - Anchor is selected among SARSA methods using global AvgRank, breaking
%   ties by GlobalMeanError, exactly as in the CEC reference logic.
% =========================================================================

clc;
try
    figs = findall(groot,'Type','figure');
    delete(figs);
catch ME
    warning('Figure cleanup skipped: %s', ME.message);
end

%% ========================= VERSIONING =========================
CFG = struct();
CFG.VERSION_TAG   = 'common_seeds';
CFG.EXP_NAME      = ['CEC2022_ablation_BC1_' CFG.VERSION_TAG];
CFG.OUTDIR        = fullfile(fileparts(mfilename('fullpath')), '..', 'results', CFG.EXP_NAME);
CFG.CLEAR_OUTDIR  = false;

if CFG.CLEAR_OUTDIR && exist(CFG.OUTDIR, 'dir')
    rmdir(CFG.OUTDIR, 's');
end
mkdir_if_needed(CFG.OUTDIR);
mkdir_if_needed(fullfile(CFG.OUTDIR, 'csv'));
mkdir_if_needed(fullfile(CFG.OUTDIR, 'tables'));
mkdir_if_needed(fullfile(CFG.OUTDIR, 'mat'));
mkdir_if_needed(fullfile(CFG.OUTDIR, 'logs'));
mkdir_if_needed(fullfile(CFG.OUTDIR, 'figures'));

CFG.SELF_COPY = true;
if CFG.SELF_COPY
    try
        thisFile = mfilename('fullpath');
        copyfile([thisFile '.m'], fullfile(CFG.OUTDIR, 'logs', [mfilename '.m']));
    catch
    end
end

%% ========================= GLOBAL SETTINGS =========================
CFG.Dim       = 20;
CFG.PopSize   = 30;
CFG.Runs      = 30;
CFG.Max_FEs   = 10000 * CFG.Dim;
CFG.Functions = 1:12;

CFG.USE_PARFOR = true;
CFG.N_WORKERS  = 8;

CFG.EXPORT_CSV   = true;
CFG.EXPORT_TEX   = true;
CFG.EXPORT_FIG   = true;
CFG.SAVE_RAW_MAT = true;

CFG.MASTER_SEED = 20260405;
CFG.TIE_TOL     = 1e-12;

CFG.RESUME_IF_AVAILABLE     = true;
CFG.RETRY_FAILED_TASKS      = true;
CFG.CHECKPOINT_BY_FUNCTION  = true;
CFG.SAVE_PARTIAL_GLOBAL     = true;
CFG.PRINT_FUNCTION_PROGRESS = true;
CFG.PROGRESS_EVERY          = 25;

CFG.FIXED_K  = 7;
CFG.FIXED_NI = 3;
CFG.T_PROBE1 = 1;
CFG.T_PROBE2 = 2;

%% ========================= ABLATION SELECTION =========================
% Best-ranked SARSA setting and its matched RANDOM mirror.
CFG.ANCHOR_CFG_ID          = 'C1';
CFG.ANCHOR_FAD             = 'BEFORE';
CFG.MATCHED_RANDOM_CFG_ID  = 'C1';
CFG.MATCHED_RANDOM_FAD     = 'BEFORE';
CFG.BEST_RANDOM_CFG_ID     = 'C2';
CFG.BEST_RANDOM_FAD        = 'BEFORE';

USER_CONFIGS = struct([]);

%% ========================= CEC2022 DEFINITION =========================
Optimos = [300, 400, 600, 800, 900, 1800, 2000, 2200, 2300, 2400, 2600, 2700];
LB = -100;
UB = 100;
cec_mex = @cec22_test_func;

%% ========================= SANITY CHECKS =========================
needFiles = {'cec22_test_func', ...
             'MPHBS_main', ...
             'MPHBS_random_mirror', ...
             'MPHB_baseline', 'HBA', 'MPA'};
for k = 1:numel(needFiles)
    if isempty(which(needFiles{k}))
        error('Required function not found on MATLAB path: %s', needFiles{k});
    end
end

%% ========================= BUILD ENTRIES =========================
Entries = build_entries_full(USER_CONFIGS, CFG);
numAlgs = numel(Entries);
numF    = numel(CFG.Functions);
numR    = CFG.Runs;

save(fullfile(CFG.OUTDIR, 'mat', 'config_snapshot.mat'), 'CFG', 'USER_CONFIGS', 'Entries', 'Optimos');
write_entries_manifest(Entries, fullfile(CFG.OUTDIR, 'csv', 'entries_manifest.csv'));

fprintf('\n============================================================\n');
fprintf('RUNNING %s\n', CFG.EXP_NAME);
fprintf('Algorithms/entries: %d\n', numAlgs);
fprintf('Functions: %d\n', numF);
fprintf('Runs per function: %d\n', numR);
fprintf('Dim: %d | PopSize: %d | Max_FEs target: %d\n', CFG.Dim, CFG.PopSize, CFG.Max_FEs);
fprintf('============================================================\n\n');

%% ========================= PARALLEL POOL =========================
if CFG.USE_PARFOR
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', CFG.N_WORKERS);
    elseif pool.NumWorkers ~= CFG.N_WORKERS
        delete(pool);
        parpool('local', CFG.N_WORKERS);
    end
end

%% ========================= GLOBAL ALLOCATION =========================
RawBest      = nan(numAlgs, numF, numR);
RawError     = nan(numAlgs, numF, numR);
RawTimeSec   = nan(numAlgs, numF, numR);
RawEstFEs    = nan(numAlgs, numF, numR);
RawCurve     = cell(numAlgs, numF, numR);
RawHitBest   = false(numAlgs, numF, numR);
RawSeedUsed  = nan(numAlgs, numF, numR);
RawOK        = false(numAlgs, numF, numR);
RawMsg       = strings(numAlgs, numF, numR);
DoneMask     = false(numAlgs, numF, numR);

%% ========================= EXECUTION BY FUNCTION =========================
for f = 1:numF
    funId = CFG.Functions(f);
    cpfile = fullfile(CFG.OUTDIR, 'mat', sprintf('checkpoint_F%d.mat', funId));

    [rawBest_f, rawErr_f, rawTime_f, rawFE_f, rawCurve_f, rawHit_f, rawSeed_f, rawOK_f, rawMsg_f, done_f] = ...
        init_or_load_function_checkpoint(cpfile, numAlgs, numR, CFG);

    if CFG.PRINT_FUNCTION_PROGRESS
        fprintf('------------------------------------------------------------\n');
        fprintf('Function F%d\n', funId);
        fprintf('Completed tasks found in checkpoint: %d / %d\n', nnz(done_f), numAlgs*numR);
    end

    missing = find(~done_f);

    if isempty(missing)
        if CFG.PRINT_FUNCTION_PROGRESS
            fprintf('F%d already complete. Loading checkpoint only.\n', funId);
        end
    else
        LocalTaskTable = build_local_task_table(numAlgs, numR, missing);
        numLocalTasks = size(LocalTaskTable, 1);
        ResultsCell = cell(numLocalTasks, 1);

        if CFG.PRINT_FUNCTION_PROGRESS
            fprintf('F%d missing tasks to execute: %d\n', funId, numLocalTasks);
        end

        if CFG.USE_PARFOR
            parfor taskId = 1:numLocalTasks
                ResultsCell{taskId} = run_one_local_task( ...
                    taskId, LocalTaskTable, Entries, CFG, LB, UB, Optimos, cec_mex, funId, f);
            end
        else
            for taskId = 1:numLocalTasks
                ResultsCell{taskId} = run_one_local_task( ...
                    taskId, LocalTaskTable, Entries, CFG, LB, UB, Optimos, cec_mex, funId, f);
            end
        end

        for taskId = 1:numLocalTasks
            rr = ResultsCell{taskId};
            a = rr.algIdx;
            r = rr.runIdx;

            rawBest_f(a,r)    = rr.bestRaw;
            rawErr_f(a,r)     = rr.err;
            rawTime_f(a,r)    = rr.timeSec;
            rawFE_f(a,r)      = rr.estFEs;
            rawCurve_f{a,r}   = rr.conv;
            rawHit_f(a,r)     = rr.hitBest;
            rawSeed_f(a,r)    = rr.seed;
            rawOK_f(a,r)      = rr.ok;
            rawMsg_f(a,r)     = string(rr.msg);
            done_f(a,r)       = rr.ok && ~isnan(rr.err);
        end

        save_function_checkpoint(cpfile, funId, f, Entries, rawBest_f, rawErr_f, rawTime_f, rawFE_f, ...
            rawCurve_f, rawHit_f, rawSeed_f, rawOK_f, rawMsg_f, done_f);

        if CFG.PRINT_FUNCTION_PROGRESS
            fprintf('F%d checkpoint saved: %d / %d tasks complete\n', funId, nnz(done_f), numAlgs*numR);
        end
    end

    RawBest(:,f,:)     = rawBest_f;
    RawError(:,f,:)    = rawErr_f;
    RawTimeSec(:,f,:)  = rawTime_f;
    RawEstFEs(:,f,:)   = rawFE_f;
    RawCurve(:,f,:)    = rawCurve_f;
    RawHitBest(:,f,:)  = rawHit_f;
    RawSeedUsed(:,f,:) = rawSeed_f;
    RawOK(:,f,:)       = rawOK_f;
    RawMsg(:,f,:)      = rawMsg_f;
    DoneMask(:,f,:)    = done_f;

    if CFG.SAVE_PARTIAL_GLOBAL
        save(fullfile(CFG.OUTDIR, 'mat', 'partial_workspace.mat'), ...
            'CFG', 'USER_CONFIGS', 'Entries', 'Optimos', 'RawBest', 'RawError', 'RawTimeSec', ...
            'RawEstFEs', 'RawCurve', 'RawHitBest', 'RawSeedUsed', 'RawOK', 'RawMsg', 'DoneMask', '-v7.3');
    end
end

%% ========================= TASK RESULTS TABLE =========================
TaskResults = rebuild_task_results_from_arrays(RawBest, RawError, RawTimeSec, RawEstFEs, RawSeedUsed, RawOK, RawMsg, Entries, CFG.Functions);
if CFG.EXPORT_CSV
    writetable(TaskResults, fullfile(CFG.OUTDIR, 'csv', 'task_results.csv'));
end

%% ========================= SUMMARY AND STATISTICS =========================
RanksRunBlock = nan(numAlgs, numF, numR);
for f = 1:numF
    for r = 1:numR
        vec = RawError(:,f,r);
        if all(isnan(vec)), continue; end
        RanksRunBlock(:,f,r) = tiedrank_local(vec')';
    end
end

Summary = build_summary_table(Entries, RawBest, RawError, RawTimeSec, RawEstFEs, RawHitBest, RawOK, DoneMask, RanksRunBlock);
[Y_valid, p_friedman, MeanRanksFromFriedman] = build_friedman_from_errors(RawError);
for a = 1:numAlgs
    Summary.FriedmanMeanRank(a) = MeanRanksFromFriedman(a);
end

[MeanByF, StdByF, RankByMeanF, AvgRankByMean] = build_functionwise_stats_from_mean(RawError);
Summary.AvgRankByMean = AvgRankByMean;
Summary = sortrows(Summary, 'FriedmanMeanRank', 'ascend');

anchorAliasTarget = build_alias('ANCHOR', CFG.ANCHOR_FAD, CFG.ANCHOR_CFG_ID);
anchorIdx = find(strcmpi(string({Entries.alias}), string(anchorAliasTarget)), 1, 'first');
if isempty(anchorIdx), error('Could not find requested anchor alias: %s', anchorAliasTarget); end
anchorName = Entries(anchorIdx).name;
anchorAlias = Entries(anchorIdx).alias;

Pairwise        = build_pairwise_vs_anchor(Entries, RawError, anchorIdx);
WTL_by_function = build_wtl_by_function(RawError, Entries, anchorIdx, CFG.Functions, CFG.TIE_TOL);
WTL_global      = build_wtl_global(RawError, Entries, anchorIdx, CFG.TIE_TOL);

%% ========================= BASE EXPORTS =========================
if CFG.EXPORT_CSV
    writetable(Summary,         fullfile(CFG.OUTDIR, 'csv', 'summary_global.csv'));
    writetable(Pairwise,        fullfile(CFG.OUTDIR, 'csv', 'pairwise_vs_anchor.csv'));
    writetable(WTL_by_function, fullfile(CFG.OUTDIR, 'csv', 'wtl_by_function_vs_anchor.csv'));
    writetable(WTL_global,      fullfile(CFG.OUTDIR, 'csv', 'wtl_global_vs_anchor.csv'));
end

%% ========================= MAIN-TEXT TABLES =========================
export_cec_ablation_tables(Entries, Pairwise, WTL_global, ...
    anchorIdx, string(anchorAlias), CFG.Functions, MeanByF, StdByF, RankByMeanF, AvgRankByMean, ...
    fullfile(CFG.OUTDIR, 'tables'), fullfile(CFG.OUTDIR, 'csv'));

%% ========================= FIGURES =========================
if CFG.EXPORT_FIG
    export_convergence_panels_cec(RawCurve, RawEstFEs, Entries, CFG, fullfile(CFG.OUTDIR, 'figures'));
end

%% ========================= SAVE FINAL MAT =========================
if CFG.SAVE_RAW_MAT
    save(fullfile(CFG.OUTDIR, 'mat', 'final_workspace.mat'), ...
        'CFG', 'USER_CONFIGS', 'Entries', 'Summary', 'Pairwise', 'WTL_by_function', ...
        'WTL_global', 'TaskResults', 'RawBest', 'RawError', 'RawTimeSec', 'RawEstFEs', ...
        'RawCurve', 'RawHitBest', 'RawSeedUsed', 'RawOK', 'RawMsg', 'DoneMask', ...
        'RanksRunBlock', 'MeanByF', 'StdByF', 'RankByMeanF', 'AvgRankByMean', ...
        'Y_valid', 'p_friedman', 'anchorIdx', 'anchorName', 'anchorAlias', '-v7.3');
end

fprintf('\n============================================================\n');
fprintf('%s DONE\n', CFG.EXP_NAME);
fprintf('Anchor for pairwise tests: %s (%s)\n', anchorAlias, anchorName);
fprintf('Friedman p-value = %.6g\n', p_friedman);
disp(Summary(:, {'Alias','Algorithm','Family','Role','GlobalMeanError','GlobalStdError','GlobalMedianError','GlobalMeanTimeSec','GlobalMeanEstFEs','AvgRankByMean','FriedmanMeanRank'}));
fprintf('Outputs saved in: %s\n', CFG.OUTDIR);
fprintf('============================================================\n');

end

%% ========================================================================
% ENTRY BUILDERS
% ========================================================================
function Entries = build_entries_full(USER_CONFIGS, CFG) %#ok<INUSD>
Entries = struct('pair_id', {}, 'cfg_id', {}, 'name', {}, 'family', {}, 'role', {}, ...
    'K', {}, 'Ni', {}, 'gamma', {}, 'w_best', {}, 'rho_sub', {}, 'label', {}, ...
    'fad_name', {}, 'fad_before_p2', {}, 'fad_after_p2', {}, 'T', {}, 'EstimatedFEs', {}, ...
    'alias', {});

idx = 0;
idx = idx + 1;
Entries(idx) = make_entry(CFG.ANCHOR_CFG_ID, 'ANCHOR_SARSA', ...
    'MPHBS_SARSA', 'ANCHOR', 'anchor_sarsa', ...
    ablation_cfg_param('K', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    CFG.ANCHOR_FAD, strcmpi(CFG.ANCHOR_FAD,'BEFORE'), strcmpi(CFG.ANCHOR_FAD,'AFTER'), ...
    build_alias('ANCHOR', CFG.ANCHOR_FAD, CFG.ANCHOR_CFG_ID));
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);

idx = idx + 1;
Entries(idx) = make_entry(CFG.MATCHED_RANDOM_CFG_ID, 'MATCHED_RANDOM', ...
    'MPHBS_RANDOM', 'MATCHED_RANDOM', 'matched_random', ...
    ablation_cfg_param('K', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    CFG.MATCHED_RANDOM_FAD, strcmpi(CFG.MATCHED_RANDOM_FAD,'BEFORE'), strcmpi(CFG.MATCHED_RANDOM_FAD,'AFTER'), ...
    build_alias('MATCHED_RANDOM', CFG.MATCHED_RANDOM_FAD, CFG.MATCHED_RANDOM_CFG_ID));
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);

idx = idx + 1;
Entries(idx) = make_entry(CFG.BEST_RANDOM_CFG_ID, 'BEST_RANDOM', ...
    'MPHBS_RANDOM', 'BEST_RANDOM', 'best_random', ...
    ablation_cfg_param('K', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    CFG.BEST_RANDOM_FAD, strcmpi(CFG.BEST_RANDOM_FAD,'BEFORE'), strcmpi(CFG.BEST_RANDOM_FAD,'AFTER'), ...
    'R-best');
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);

idx = idx + 1;
Entries(idx) = make_entry('BASE', 'BASE', 'MPHB_BASE', 'BASELINE', 'baseline', ...
    nan, nan, nan, nan, nan, 'BASE', false, false, 'MPHB');
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);

idx = idx + 1;
Entries(idx) = make_entry('HBA', 'HBA', 'HBA_BASE', 'PARENT', 'parent_hba', ...
    nan, nan, nan, nan, nan, 'BASE', false, false, 'HBA');
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);

idx = idx + 1;
Entries(idx) = make_entry('MPA', 'MPA', 'MPA_BASE', 'PARENT', 'parent_mpa', ...
    nan, nan, nan, nan, nan, 'BASE', false, false, 'MPA');
[Entries(idx).T, Entries(idx).EstimatedFEs] = estimate_T_and_FEs(Entries(idx), CFG);
end

function v = ablation_cfg_param(whichParam, cfg_id, fad_name, CFG) %#ok<INUSD>
switch upper(whichParam)
    case 'K'
        v = CFG.FIXED_K;
    case 'NI'
        v = CFG.FIXED_NI;
    case 'W_BEST'
        switch upper(char(cfg_id))
            case 'C1', v = 0.050;
            case 'C2', v = 0.025;
            case 'C3', v = 0.050;
            case 'C4', v = 0.050;
            otherwise, error('Unknown cfg id: %s', char(cfg_id));
        end
    case 'GAMMA'
        switch upper(char(cfg_id))
            case 'C1', v = 0.30;
            case 'C2', v = 0.30;
            case 'C3', v = 0.20;
            case 'C4', v = 0.30;
            otherwise, error('Unknown cfg id: %s', char(cfg_id));
        end
    case 'RHO_SUB'
        switch upper(char(cfg_id))
            case {'C1','C2','C3'}, v = 0.25;
            case 'C4', v = 0.30;
            otherwise, error('Unknown cfg id: %s', char(cfg_id));
        end
    otherwise
        error('Unknown parameter request: %s', whichParam);
end
end

function cfg = make_cfg(cfg_id, w_best, gamma, rho_sub, fad_name, enabled) %#ok<INUSD>
cfg = struct();
cfg.cfg_id         = string(cfg_id);
cfg.label          = string(cfg_id) + "_" + string(fad_name);
cfg.w_best         = w_best;
cfg.gamma          = gamma;
cfg.rho_sub        = rho_sub;
cfg.fad_name       = string(fad_name);
cfg.fad_before_p2  = strcmpi(fad_name, 'BEFORE');
cfg.fad_after_p2   = strcmpi(fad_name, 'AFTER');
cfg.enabled        = logical(enabled);
end

function e = make_entry(cfg_id, pair_id, family, role, label, K, Ni, gamma, w_best, rho_sub, fad_name, fad_before_p2, fad_after_p2, alias)
e = struct();
e.cfg_id         = char(cfg_id);
e.pair_id        = char(pair_id);
e.family         = char(family);
e.role           = char(role);
e.label          = char(label);
e.K              = K;
e.Ni             = Ni;
e.gamma          = gamma;
e.w_best         = w_best;
e.rho_sub        = rho_sub;
e.fad_name       = char(fad_name);
e.fad_before_p2  = logical(fad_before_p2);
e.fad_after_p2   = logical(fad_after_p2);
e.T              = nan;
e.EstimatedFEs   = nan;
e.alias          = char(alias);
e.name           = build_entry_name(e);
end

function alias = build_alias(role, fad_name, cfg_id)
cfgnum = replace(string(cfg_id), "C", "");
switch upper(string(role))
    case "ANCHOR"
        if strcmpi(fad_name, 'BEFORE')
            alias = "BS-C" + cfgnum;
        else
            alias = "AS-C" + cfgnum;
        end
    case "MATCHED_RANDOM"
        if strcmpi(fad_name, 'BEFORE')
            alias = "BR-C" + cfgnum;
        else
            alias = "AR-C" + cfgnum;
        end
    otherwise
        if strcmpi(fad_name, 'BEFORE')
            alias = "B-C" + cfgnum;
        else
            alias = "A-C" + cfgnum;
        end
end
end

function name = build_entry_name(e)
switch upper(e.family)
    case 'MPHBS_SARSA'
        name = sprintf('MPHBS_SARSA_%s_K%02dN%02d_%s_W%03d_G%03d_R%03d', ...
            e.cfg_id, e.K, e.Ni, upper(e.fad_name), ...
            round(1000*e.w_best), round(1000*e.gamma), round(1000*e.rho_sub));
    case 'MPHBS_RANDOM'
        name = sprintf('MPHBS_RANDOM_%s_K%02dN%02d_%s_W%03d_G%03d_R%03d', ...
            e.cfg_id, e.K, e.Ni, upper(e.fad_name), ...
            round(1000*e.w_best), round(1000*e.gamma), round(1000*e.rho_sub));
    case 'MPHB_BASE'
        name = 'MPHB_BASE';
    case 'HBA_BASE'
        name = 'HBA_BASE';
    case 'MPA_BASE'
        name = 'MPA_BASE';
    otherwise
        name = char(e.family);
end
end

function [T, estFE] = estimate_T_and_FEs(entry, CFG)
probeDim = CFG.Dim;
T1 = CFG.T_PROBE1;
T2 = CFG.T_PROBE2;
fe1 = probe_actual_fes(entry, CFG, probeDim, T1);
fe2 = probe_actual_fes(entry, CFG, probeDim, T2);
slope = fe2 - fe1;
intercept = fe1 - slope * T1;
if ~(isfinite(slope) && isfinite(intercept)) || slope <= 0
    error('Invalid FE probe for family=%s, cfg=%s (fe1=%g, fe2=%g).', entry.family, entry.cfg_id, fe1, fe2);
end
T = floor((CFG.Max_FEs - intercept) / slope);
T = max(1, T);
estFE = intercept + slope * T;
while estFE > CFG.Max_FEs && T > 1
    T = T - 1;
    estFE = intercept + slope * T;
end
end

function fe_used = probe_actual_fes(entry, CFG, probeDim, Tprobe)
global FE_COUNTER_CEC_ABL
FE_COUNTER_CEC_ABL = 0;
rng(12345, 'twister');
f_raw = @(x) sum(x(:).^2);
F = @(x) probe_count_objective_cec(x, f_raw);
switch upper(entry.family)
    case 'MPHBS_SARSA'
        MPHBS_main(CFG.PopSize, Tprobe, -100, 100, probeDim, F, entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, entry.fad_before_p2, entry.fad_after_p2);
    case 'MPHBS_RANDOM'
        MPHBS_random_mirror(CFG.PopSize, Tprobe, -100, 100, probeDim, F, entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, entry.fad_before_p2, entry.fad_after_p2);
    case 'MPHB_BASE'
        MPHB_baseline(CFG.PopSize, Tprobe, -100, 100, probeDim, F);
    case 'HBA_BASE'
        HBA(F, probeDim, -100, 100, Tprobe, CFG.PopSize);
    case 'MPA_BASE'
        MPA(CFG.PopSize, Tprobe, -100, 100, probeDim, F);
    otherwise
        error('Unknown family for FE probe: %s', entry.family);
end
fe_used = FE_COUNTER_CEC_ABL;
end

function y = probe_count_objective_cec(x, f_raw)
global FE_COUNTER_CEC_ABL
x = x(:).';
FE_COUNTER_CEC_ABL = FE_COUNTER_CEC_ABL + 1;
y = f_raw(x);
end

function write_entries_manifest(Entries, outcsv)
n = numel(Entries);
T = table();
T.Alias     = string({Entries.alias})';
T.Algorithm = string({Entries.name})';
T.Family    = string({Entries.family})';
T.Role      = string({Entries.role})';
T.CFG_ID    = string({Entries.cfg_id})';
T.FAD_Mode  = string({Entries.fad_name})';
Kvec = nan(n,1); Nivec = nan(n,1); Gvec = nan(n,1); Wvec = nan(n,1); Rvec = nan(n,1); Tvec = nan(n,1); FEvec = nan(n,1);
for i = 1:n
    Kvec(i)=Entries(i).K; Nivec(i)=Entries(i).Ni; Gvec(i)=Entries(i).gamma;
    Wvec(i)=Entries(i).w_best; Rvec(i)=Entries(i).rho_sub; Tvec(i)=Entries(i).T; FEvec(i)=Entries(i).EstimatedFEs;
end
T.K = Kvec; T.Ni = Nivec; T.gamma = Gvec; T.w_best = Wvec; T.rho_sub = Rvec; T.T = Tvec; T.EstimatedFEs = FEvec;
writetable(T, outcsv);
end

%% ========================================================================
% EXECUTION HELPERS
% ========================================================================
function rr = run_one_local_task(taskId, LocalTaskTable, Entries, CFG, LB, UB, Optimos, cec_mex, funId, function_pos)
a = LocalTaskTable(taskId,1);
r = LocalTaskTable(taskId,2);
entry = Entries(a);

% Common random numbers across ablation methods for this function/run block.
seed = CFG.MASTER_SEED + 1000*function_pos + r;
rng(seed, 'twister');

global FE_COUNTER_CEC_REAL
FE_COUNTER_CEC_REAL = 0;
fobj = @(x) cec_error_wrapper(x, funId, cec_mex, Optimos);

try
    t0 = tic;
    switch upper(entry.family)
        case 'MPHBS_SARSA'
            [Pg, ~, Conv] = MPHBS_main( ...
                CFG.PopSize, entry.T, LB, UB, CFG.Dim, fobj, ...
                entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, ...
                entry.fad_before_p2, entry.fad_after_p2);

        case 'MPHBS_RANDOM'
            [Pg, ~, Conv] = MPHBS_random_mirror( ...
                CFG.PopSize, entry.T, LB, UB, CFG.Dim, fobj, ...
                entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, ...
                entry.fad_before_p2, entry.fad_after_p2);

        case 'MPHB_BASE'
            [Pg, ~, Conv] = MPHB_baseline( ...
                CFG.PopSize, entry.T, LB, UB, CFG.Dim, fobj);

        case 'HBA_BASE'
            [bestX, Pg, Conv] = HBA(fobj, CFG.Dim, LB, UB, entry.T, CFG.PopSize); %#ok<ASGLU>

        case 'MPA_BASE'
            [Pg, bestX, Conv] = MPA(CFG.PopSize, entry.T, LB, UB, CFG.Dim, fobj); %#ok<ASGLU>

        otherwise
            error('Unknown family: %s', entry.family);
    end
    elapsed = toc(t0);

    rr.algIdx  = a;
    rr.runIdx  = r;
    rr.seed    = seed;
    rr.bestRaw = Pg + Optimos(funId);
    rr.err     = Pg;
    rr.timeSec = elapsed;
    rr.estFEs  = FE_COUNTER_CEC_REAL;
    rr.conv    = force_row_curve(Conv);
    rr.hitBest = abs(Pg) <= 1e-8;
    rr.ok      = true;
    rr.msg     = "OK";
catch ME
    rr.algIdx  = a;
    rr.runIdx  = r;
    rr.seed    = seed;
    rr.bestRaw = nan;
    rr.err     = nan;
    rr.timeSec = nan;
    rr.estFEs  = FE_COUNTER_CEC_REAL;
    rr.conv    = [];
    rr.hitBest = false;
    rr.ok      = false;
    rr.msg     = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
end
end

function err = cec_error_wrapper(x, funId, cec_mex, Optimos)
global FE_COUNTER_CEC_REAL
FE_COUNTER_CEC_REAL = FE_COUNTER_CEC_REAL + 1;
x = x(:)';
ok = false;
try
    fraw = cec_mex(x', funId);
    ok = true;
catch
end
if ~ok
    fraw = cec_mex(x, funId);
end
if numel(fraw) > 1
    fraw = fraw(1);
end
err = fraw - Optimos(funId);
end

function y = force_row_curve(x)
if isempty(x), y = []; else, y = x(:).'; end
end

function [rawBest_f, rawErr_f, rawTime_f, rawFE_f, rawCurve_f, rawHit_f, rawSeed_f, rawOK_f, rawMsg_f, done_f] = ...
    init_or_load_function_checkpoint(cpfile, numAlgs, numR, CFG)

rawBest_f = nan(numAlgs, numR);
rawErr_f  = nan(numAlgs, numR);
rawTime_f = nan(numAlgs, numR);
rawFE_f   = nan(numAlgs, numR);
rawCurve_f = cell(numAlgs, numR);
rawHit_f  = false(numAlgs, numR);
rawSeed_f = nan(numAlgs, numR);
rawOK_f   = false(numAlgs, numR);
rawMsg_f  = strings(numAlgs, numR);
done_f    = false(numAlgs, numR);

if ~CFG.RESUME_IF_AVAILABLE || ~exist(cpfile, 'file')
    return;
end

try
    S = load(cpfile);
    if isfield(S, 'rawBest_f') && isequal(size(S.rawBest_f), [numAlgs, numR]), rawBest_f = S.rawBest_f; end
    if isfield(S, 'rawErr_f')  && isequal(size(S.rawErr_f),  [numAlgs, numR]), rawErr_f  = S.rawErr_f;  end
    if isfield(S, 'rawTime_f') && isequal(size(S.rawTime_f), [numAlgs, numR]), rawTime_f = S.rawTime_f; end
    if isfield(S, 'rawFE_f')   && isequal(size(S.rawFE_f),   [numAlgs, numR]), rawFE_f   = S.rawFE_f;   end
    if isfield(S, 'rawCurve_f')&& isequal(size(S.rawCurve_f),[numAlgs, numR]), rawCurve_f= S.rawCurve_f;end
    if isfield(S, 'rawHit_f')  && isequal(size(S.rawHit_f),  [numAlgs, numR]), rawHit_f  = S.rawHit_f;  end
    if isfield(S, 'rawSeed_f') && isequal(size(S.rawSeed_f), [numAlgs, numR]), rawSeed_f = S.rawSeed_f; end
    if isfield(S, 'rawOK_f')   && isequal(size(S.rawOK_f),   [numAlgs, numR]), rawOK_f   = S.rawOK_f;   end
    if isfield(S, 'rawMsg_f')  && isequal(size(S.rawMsg_f),  [numAlgs, numR]), rawMsg_f  = string(S.rawMsg_f); end
    if isfield(S, 'done_f') && isequal(size(S.done_f), [numAlgs, numR])
        done_f = S.done_f;
    else
        done_f = ~isnan(rawErr_f);
    end

    if CFG.RETRY_FAILED_TASKS
        bad = (~rawOK_f) | isnan(rawErr_f);
        done_f(bad) = false;
        rawBest_f(bad) = nan;
        rawErr_f(bad)  = nan;
        rawTime_f(bad) = nan;
        rawFE_f(bad)   = nan;
        rawHit_f(bad)  = false;
        rawSeed_f(bad) = nan;
        rawMsg_f(bad)  = "";
        for ii = 1:numel(bad)
            if bad(ii), rawCurve_f{ii} = []; end
        end
    end
catch
end
end

function save_function_checkpoint(cpfile, funId, funcPos, Entries, rawBest_f, rawErr_f, rawTime_f, rawFE_f, rawCurve_f, rawHit_f, rawSeed_f, rawOK_f, rawMsg_f, done_f)
save(cpfile, 'funId', 'funcPos', 'Entries', 'rawBest_f', 'rawErr_f', 'rawTime_f', 'rawFE_f', ...
    'rawCurve_f', 'rawHit_f', 'rawSeed_f', 'rawOK_f', 'rawMsg_f', 'done_f', '-v7.3');
end

function LocalTaskTable = build_local_task_table(numAlgs, numR, linIdxMissing)
LocalTaskTable = zeros(numel(linIdxMissing), 3);
for k = 1:numel(linIdxMissing)
    [a, r] = ind2sub([numAlgs, numR], linIdxMissing(k));
    LocalTaskTable(k,:) = [a, r, linIdxMissing(k)];
end
end

function TaskResults = rebuild_task_results_from_arrays(RawBest, RawError, RawTimeSec, RawEstFEs, RawSeedUsed, RawOK, RawMsg, Entries, function_list)
numAlgs = size(RawError,1);
numF    = size(RawError,2);
numR    = size(RawError,3);
numTasks = numAlgs * numF * numR;

TaskResults = table('Size', [numTasks 16], ...
    'VariableTypes', {'double','double','double','double','string','string','string','string','string','double','logical','double','double','double','double','string'}, ...
    'VariableNames', {'AlgIdx','FuncPos','FunctionID','RunIdx','Alias','Algorithm','Family','Role','PairID','Seed','OK','BestRaw','Error','TimeSec','EstFEs','Message'});

row = 0;
for a = 1:numAlgs
    for f = 1:numF
        for r = 1:numR
            row = row + 1;
            TaskResults.AlgIdx(row)     = a;
            TaskResults.FuncPos(row)    = f;
            TaskResults.FunctionID(row) = function_list(f);
            TaskResults.RunIdx(row)     = r;
            TaskResults.Alias(row)      = string(Entries(a).alias);
            TaskResults.Algorithm(row)  = string(Entries(a).name);
            TaskResults.Family(row)     = string(Entries(a).family);
            TaskResults.Role(row)       = string(Entries(a).role);
            TaskResults.PairID(row)     = string(Entries(a).pair_id);
            TaskResults.Seed(row)       = RawSeedUsed(a,f,r);
            TaskResults.OK(row)         = RawOK(a,f,r);
            TaskResults.BestRaw(row)    = RawBest(a,f,r);
            TaskResults.Error(row)      = RawError(a,f,r);
            TaskResults.TimeSec(row)    = RawTimeSec(a,f,r);
            TaskResults.EstFEs(row)     = RawEstFEs(a,f,r);
            TaskResults.Message(row)    = string(RawMsg(a,f,r));
        end
    end
end
end

%% ========================================================================
% SUMMARY / STATS
% ========================================================================
function Summary = build_summary_table(Entries, RawBest, RawError, RawTimeSec, RawEstFEs, RawHitBest, RawOK, DoneMask, RanksRunBlock)
numAlgs = numel(Entries);
Summary = table();
Summary.Alias             = strings(numAlgs,1);
Summary.Algorithm         = strings(numAlgs,1);
Summary.Family            = strings(numAlgs,1);
Summary.Role              = strings(numAlgs,1);
Summary.GlobalMeanRaw     = nan(numAlgs,1);
Summary.GlobalMeanError   = nan(numAlgs,1);
Summary.GlobalStdError    = nan(numAlgs,1);
Summary.GlobalMedianError = nan(numAlgs,1);
Summary.GlobalMeanTimeSec = nan(numAlgs,1);
Summary.GlobalMeanEstFEs  = nan(numAlgs,1);
Summary.SuccessRatePct    = nan(numAlgs,1);
Summary.ValidPct          = nan(numAlgs,1);
Summary.CompletedPct      = nan(numAlgs,1);
Summary.MeanRankRunFunc   = nan(numAlgs,1);
Summary.FriedmanMeanRank  = nan(numAlgs,1);

for a = 1:numAlgs
    Summary.Alias(a)             = string(Entries(a).alias);
    Summary.Algorithm(a)         = string(Entries(a).name);
    Summary.Family(a)            = string(Entries(a).family);
    Summary.Role(a)              = string(Entries(a).role);
    Summary.GlobalMeanRaw(a)     = mean(reshape(RawBest(a,:,:), 1, []), 'omitnan');
    Summary.GlobalMeanError(a)   = mean(reshape(RawError(a,:,:),1, []), 'omitnan');
    Summary.GlobalStdError(a)    = std(reshape(RawError(a,:,:),1, []), 0, 'omitnan');
    Summary.GlobalMedianError(a) = median(reshape(RawError(a,:,:),1, []), 'omitnan');
    Summary.GlobalMeanTimeSec(a) = mean(reshape(RawTimeSec(a,:,:),1, []), 'omitnan');
    Summary.GlobalMeanEstFEs(a)  = mean(reshape(RawEstFEs(a,:,:),1, []), 'omitnan');
    Summary.SuccessRatePct(a)    = 100 * mean(reshape(RawHitBest(a,:,:),1, []), 'omitnan');
    Summary.ValidPct(a)          = 100 * mean(reshape(RawOK(a,:,:),1, []), 'omitnan');
    Summary.CompletedPct(a)      = 100 * mean(reshape(DoneMask(a,:,:),1, []), 'omitnan');
    Summary.MeanRankRunFunc(a)   = mean(reshape(RanksRunBlock(a,:,:),1, []), 'omitnan');
end
end

function [Y_valid, p_friedman, MeanRanksFromFriedman] = build_friedman_from_errors(RawError)
[numAlgs, numF, numR] = size(RawError);
Y = nan(numF*numR, numAlgs);
row = 0;
for f = 1:numF
    for r = 1:numR
        row = row + 1;
        Y(row,:) = RawError(:,f,r)';
    end
end
validBlocks = all(~isnan(Y), 2);
Y_valid = Y(validBlocks, :);

if size(Y_valid,1) >= 2 && size(Y_valid,2) >= 2
    [p_friedman, ~, friedmanStats] = friedman(Y_valid, 1, 'off');
    MeanRanksFromFriedman = friedmanStats.meanranks(:);
else
    p_friedman = nan;
    MeanRanksFromFriedman = nan(numAlgs,1);
end
end

function [MeanByF, StdByF, RankByMeanF, AvgRankByMean] = build_functionwise_stats_from_mean(RawError)
[numAlgs, numF, ~] = size(RawError);
MeanByF = nan(numF, numAlgs);
StdByF  = nan(numF, numAlgs);
RankByMeanF = nan(numF, numAlgs);

for f = 1:numF
    for a = 1:numAlgs
        x = squeeze(RawError(a,f,:));
        MeanByF(f,a) = mean(x, 'omitnan');
        StdByF(f,a)  = std(x, 0, 'omitnan');
    end
    RankByMeanF(f,:) = tiedrank_local(MeanByF(f,:)')';
end

AvgRankByMean = mean(RankByMeanF, 1, 'omitnan')';
end

function anchorIdx = pick_anchor_index_sarsa_by_meanrank(Entries, AvgRankByMean, Summary)
roles = string({Entries.role})';
cand = find(roles == "SARSA");
[~, ord] = sort(AvgRankByMean(cand), 'ascend');
anchorIdx = cand(ord(1));

bestVal = AvgRankByMean(anchorIdx);
tie = cand(abs(AvgRankByMean(cand) - bestVal) <= 1e-12);
if numel(tie) > 1
    sumNames = string(Summary.Algorithm);
    vals = nan(numel(tie),1);
    for i = 1:numel(tie)
        idx = find(sumNames == string(Entries(tie(i)).name), 1, 'first');
        vals(i) = Summary.GlobalMeanError(idx);
    end
    [~, j] = min(vals);
    anchorIdx = tie(j);
end
end

function Pairwise = build_pairwise_vs_anchor(Entries, RawError, anchorIdx)
numAlgs = numel(Entries);
Pairwise = table();
Pairwise.Alias                = strings(numAlgs-1,1);
Pairwise.Algorithm            = strings(numAlgs-1,1);
Pairwise.Family               = strings(numAlgs-1,1);
Pairwise.Role                 = strings(numAlgs-1,1);
Pairwise.Anchor               = strings(numAlgs-1,1);
Pairwise.RawP                 = nan(numAlgs-1,1);
Pairwise.HolmAdjP             = nan(numAlgs-1,1);
Pairwise.RejectHolm           = false(numAlgs-1,1);
Pairwise.MedianDelta          = nan(numAlgs-1,1);
Pairwise.BetterThanAnchorPct  = nan(numAlgs-1,1);

allAnchor = reshape(RawError(anchorIdx,:,:), 1, []);
keep = 0; rawp = []; rowmap = [];
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    keep = keep + 1;
    x = reshape(RawError(a,:,:), 1, []);
    [p, ~, ~] = signrank_safe(x, allAnchor);

    Pairwise.Alias(keep)     = string(Entries(a).alias);
    Pairwise.Algorithm(keep) = string(Entries(a).name);
    Pairwise.Family(keep)    = string(Entries(a).family);
    Pairwise.Role(keep)      = string(Entries(a).role);
    Pairwise.Anchor(keep)    = string(Entries(anchorIdx).name);
    Pairwise.RawP(keep)      = p;
    Pairwise.MedianDelta(keep) = median(x - allAnchor, 'omitnan');
    Pairwise.BetterThanAnchorPct(keep) = 100 * mean(x < allAnchor, 'omitnan');

    rawp(end+1,1) = p; %#ok<AGROW>
    rowmap(end+1,1) = keep; %#ok<AGROW>
end
Pairwise = Pairwise(1:keep,:);

validP = ~isnan(rawp);
adjp_all = nan(size(rawp));
rej_all  = false(size(rawp));
if any(validP)
    [adjp_tmp, rej_tmp] = holm_correction(rawp(validP), 0.05);
    adjp_all(validP) = adjp_tmp;
    rej_all(validP)  = rej_tmp;
end
for i = 1:numel(rowmap)
    Pairwise.HolmAdjP(rowmap(i))   = adjp_all(i);
    Pairwise.RejectHolm(rowmap(i)) = rej_all(i);
end
Pairwise = sortrows(Pairwise, 'RawP', 'ascend');
end

function Tbl = build_wtl_by_function(RawError, Entries, anchorIdx, function_list, tieTol)
numAlgs = size(RawError,1);
numF    = size(RawError,2);
maxRows = (numAlgs - 1) * numF;
Tbl = table('Size', [maxRows 7], ...
    'VariableTypes', {'string','string','string','string','double','double','string'}, ...
    'VariableNames', {'Alias','Algorithm','Family','Role','Function','AnchorIndex','WTL'});
row = 0;
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    for f = 1:numF
        row = row + 1;
        x = squeeze(RawError(a,f,:));
        y = squeeze(RawError(anchorIdx,f,:));
        [w, t, l] = count_wtl(x, y, tieTol);
        Tbl.Alias(row)       = string(Entries(a).alias);
        Tbl.Algorithm(row)   = string(Entries(a).name);
        Tbl.Family(row)      = string(Entries(a).family);
        Tbl.Role(row)        = string(Entries(a).role);
        Tbl.Function(row)    = function_list(f);
        Tbl.AnchorIndex(row) = anchorIdx;
        Tbl.WTL(row)         = sprintf('%d/%d/%d', w, t, l);
    end
end
Tbl = Tbl(1:row,:);
end

function Tbl = build_wtl_global(RawError, Entries, anchorIdx, tieTol)
numAlgs = size(RawError,1);
rows = numAlgs - 1;
Tbl = table('Size', [rows 6], ...
    'VariableTypes', {'string','string','string','string','string','double'}, ...
    'VariableNames', {'Alias','Algorithm','Family','Role','WTL_Global','AnchorIndex'});
row = 0;
anchor = reshape(RawError(anchorIdx,:,:), 1, []);
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    row = row + 1;
    x = reshape(RawError(a,:,:),1,[]);
    [w, t, l] = count_wtl(x, anchor, tieTol);
    Tbl.Alias(row)       = string(Entries(a).alias);
    Tbl.Algorithm(row)   = string(Entries(a).name);
    Tbl.Family(row)      = string(Entries(a).family);
    Tbl.Role(row)        = string(Entries(a).role);
    Tbl.WTL_Global(row)  = sprintf('%d/%d/%d', w, t, l);
    Tbl.AnchorIndex(row) = anchorIdx;
end
Tbl = Tbl(1:row,:);
end

function [w, t, l] = count_wtl(x, y, tieTol)
mask = ~(isnan(x) | isnan(y));
x = x(mask); y = y(mask);
if isempty(x), w = 0; t = 0; l = 0; return; end
d = x - y;
w = sum(d < -tieTol);
t = sum(abs(d) <= tieTol);
l = sum(d > tieTol);
end

%% ========================================================================
% JOURNAL TABLE EXPORT (UNIFIED)
% ========================================================================
function export_cec_ablation_tables(Entries, Pairwise, WTL_global, ...
    anchorIdx, anchorAlias, function_list, MeanByF, StdByF, RankByMeanF, AvgRankAll, tables_dir, csv_dir)

idxKeep = 1:numel(Entries);
ItemLabels = "F" + string(function_list(:));
AliasTbl = build_alias_table_ablation(Entries);
MainTbl = build_sensitivity_panel_table_unified(Entries, idxKeep, ItemLabels, 'Function', MeanByF, StdByF, RankByMeanF, AvgRankAll);
StatsTbl = build_anchor_stats_table_ablation(Entries, Pairwise, WTL_global, AvgRankAll, anchorIdx);

writetable(AliasTbl, fullfile(csv_dir, 'cec_ablation_alias_map.csv'));
writetable(MainTbl,  fullfile(csv_dir, 'cec_ablation_main.csv'));
writetable(StatsTbl, fullfile(csv_dir, 'cec_ablation_anchor_stats.csv'));

write_latex_alias_map_ablation(AliasTbl, fullfile(tables_dir, 'cec_ablation_alias_map.tex'));
write_latex_ablation_panel_unified(Entries, idxKeep, ItemLabels, 'Function', 'CEC2022', 'error', MeanByF, StdByF, RankByMeanF, AvgRankAll, ...
    fullfile(tables_dir, 'cec_ablation_main.tex'));
write_latex_anchor_stats_unified(StatsTbl, string(anchorAlias), fullfile(tables_dir, 'cec_ablation_anchor_stats.tex'));
end

function Tbl = build_alias_table_ablation(Entries)
n = numel(Entries);
Alias = strings(n,1); Method = strings(n,1); Definition = strings(n,1);
for i = 1:n
    Alias(i) = string(Entries(i).alias);
    Method(i) = string(Entries(i).name);
    switch upper(string(Entries(i).family))
        case 'MPHBS_SARSA'
            Definition(i) = sprintf('SARSA anchor | K=%d, Ni=%d, %s | w_best=%.3f, gamma=%.2f, rho=%.2f', Entries(i).K, Entries(i).Ni, upper(Entries(i).fad_name), Entries(i).w_best, Entries(i).gamma, Entries(i).rho_sub);
        case 'MPHBS_RANDOM'
            if strcmpi(Entries(i).role,'MATCHED_RANDOM')
                prefix = 'RANDOM matched';
            else
                prefix = 'RANDOM best';
            end
            Definition(i) = sprintf('%s | K=%d, Ni=%d, %s | w_best=%.3f, gamma=%.2f, rho=%.2f', prefix, Entries(i).K, Entries(i).Ni, upper(Entries(i).fad_name), Entries(i).w_best, Entries(i).gamma, Entries(i).rho_sub);
        case 'MPHB_BASE'
            Definition(i) = 'MPHB baseline';
        case 'HBA_BASE'
            Definition(i) = 'HBA parent algorithm';
        case 'MPA_BASE'
            Definition(i) = 'MPA parent algorithm';
        otherwise
            Definition(i) = string(Entries(i).family);
    end
end
Tbl = table(Alias, Method, Definition);
end

function Tbl = build_anchor_stats_table_ablation(Entries, Pairwise, WTL_global, AvgRankAll, anchorIdx)
n = numel(Entries);
aliases = string({Entries.alias})';
algnames = string({Entries.name})';
roles = string({Entries.role})';
pairNames = string(Pairwise.Algorithm);
wtlNames = string(WTL_global.Algorithm);
idxKeep = setdiff((1:n)', anchorIdx);
Alias = strings(numel(idxKeep),1); Algorithm = strings(numel(idxKeep),1); Role = strings(numel(idxKeep),1);
GlobalWTL = strings(numel(idxKeep),1); AvgRank = nan(numel(idxKeep),1); RawP = nan(numel(idxKeep),1); HolmP = nan(numel(idxKeep),1); RejectHolm = false(numel(idxKeep),1); BetterPct = nan(numel(idxKeep),1);
for j = 1:numel(idxKeep)
    i = idxKeep(j);
    Alias(j) = aliases(i); Algorithm(j) = algnames(i); Role(j) = roles(i); AvgRank(j) = AvgRankAll(i);
    pidx = find(pairNames == algnames(i), 1, 'first');
    if ~isempty(pidx)
        RawP(j) = Pairwise.RawP(pidx); HolmP(j) = Pairwise.HolmAdjP(pidx); RejectHolm(j) = Pairwise.RejectHolm(pidx); BetterPct(j) = Pairwise.BetterThanAnchorPct(pidx);
    end
    gidx = find(wtlNames == algnames(i), 1, 'first');
    if ~isempty(gidx), GlobalWTL(j) = string(WTL_global.WTL_Global(gidx)); end
end
Tbl = table(Alias, Algorithm, Role, GlobalWTL, AvgRank, RawP, HolmP, RejectHolm, BetterPct);
Tbl = sortrows(Tbl, {'AvgRank','RawP'}, {'ascend','ascend'});
end

function write_latex_alias_map_ablation(T, fname)
fid = fopen(fname,'w'); assert(fid>0,'Could not create %s', fname);
fprintf(fid, '%% Auto-generated\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
fprintf(fid, '\\caption{Alias map for the CEC2022 ablation study.}\n');
fprintf(fid, '\\begin{tabular}{l l l}\\hline\n');
fprintf(fid, 'Alias & Method & Definition \\\\ \\hline\n');
for i = 1:height(T)
    fprintf(fid, '%s & %s & %s \\\\ \n', ...
        escape_tex(T.Alias(i)), escape_tex(T.Method(i)), escape_tex(T.Definition(i)));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
end

function write_latex_ablation_panel_unified(Entries, idxKeep, ItemLabels, itemName, domainName, metricName, MeanMat, StdMat, RankMat, AvgRankAll, fname)
fid = fopen(fname,'w'); assert(fid>0,'Could not create %s', fname);
fprintf(fid, '%% Auto-generated\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
cap = sprintf('%s ablation results for the selected anchor, its matched RANDOM mirror, the best RANDOM configuration, MPHB, HBA, and MPA. For each %s, the %s mean, standard deviation, and global rank are reported in separate rows. Global ranks are computed from the full ablation set using mean values. Lower values are better. The minimum mean within each displayed row block is highlighted in bold.', domainName, lower(itemName), metricName);
fprintf(fid, '\\caption{%s}\n', escape_tex(cap));
fprintf(fid, '\\resizebox{\\textwidth}{!}{\n');
fprintf(fid, '\\begin{tabular}{ll');
for k = 1:numel(idxKeep)
    fprintf(fid, 'c');
end
fprintf(fid, '}\\n');
fprintf(fid, '\\hline\n');
fprintf(fid, '%s & Statistic', itemName);
for k = 1:numel(idxKeep)
    fprintf(fid, ' & %s', escape_tex(string(Entries(idxKeep(k)).alias)));
end
fprintf(fid, ' \\\\ \\hline\n');
for i = 1:numel(ItemLabels)
    vals = nan(1,numel(idxKeep));
    for k = 1:numel(idxKeep)
        vals(k) = MeanMat(i, idxKeep(k));
    end
    [~, bestk] = min(vals);
    fprintf(fid, '\\multirow{3}{*}{%s} & mean', escape_tex(ItemLabels(i)));
    for k = 1:numel(idxKeep)
        if k == bestk
            fprintf(fid, ' & \\textbf{%s}', sprintf('%.3e', MeanMat(i, idxKeep(k))));
        else
            fprintf(fid, ' & %s', sprintf('%.3e', MeanMat(i, idxKeep(k))));
        end
    end
    fprintf(fid, ' \\\\ \n');
    fprintf(fid, ' & std');
    for k = 1:numel(idxKeep)
        fprintf(fid, ' & %s', sprintf('%.3e', StdMat(i, idxKeep(k))));
    end
    fprintf(fid, ' \\\\ \n');
    fprintf(fid, ' & rank');
    for k = 1:numel(idxKeep)
        fprintf(fid, ' & %s', sprintf('%.3f', RankMat(i, idxKeep(k))));
    end
    fprintf(fid, ' \\\\ \\hline\n');
end
avg = nan(1,numel(idxKeep));
for k = 1:numel(idxKeep)
    avg(k) = AvgRankAll(idxKeep(k));
end
[~, bestk] = min(avg);
fprintf(fid, 'AvgRank & rank');
for k = 1:numel(idxKeep)
    if k == bestk
        fprintf(fid, ' & \\textbf{%s}', sprintf('%.3f', AvgRankAll(idxKeep(k))));
    else
        fprintf(fid, ' & %s', sprintf('%.3f', AvgRankAll(idxKeep(k))));
    end
end
fprintf(fid, ' \\\\ \n');
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
end

function Tbl = build_alias_table_unified(Entries)
n = numel(Entries);
Alias = strings(n,1);
Method = strings(n,1);
Definition = strings(n,1);
for i = 1:n
    Alias(i) = string(Entries(i).alias);
    Method(i) = string(Entries(i).name);
    if Alias(i) == "MPHB"
        Definition(i) = "MPHB baseline";
    else
        Definition(i) = sprintf('%s | K=%d, Ni=%d, %s | w_best=%.3f, gamma=%.2f, rho=%.2f', ...
            Entries(i).role, Entries(i).K, Entries(i).Ni, upper(Entries(i).fad_name), Entries(i).w_best, Entries(i).gamma, Entries(i).rho_sub);
    end
end
Tbl = table(Alias, Method, Definition);
end

function idxKeep = build_panel_index(Entries, panelRole)
idxPanel = find(strcmpi({Entries.role}, panelRole));
idxBase  = find(strcmpi({Entries.alias}, 'MPHB'), 1, 'first');
idxKeep  = [idxPanel(:); idxBase];
end

function Tbl = build_sensitivity_panel_table_unified(Entries, idxKeep, ItemLabels, itemName, MeanMat, StdMat, RankMat, AvgRankAll)
nItems = numel(ItemLabels);
nAlg = numel(idxKeep);
nRows = 3*nItems + 1;

Item = strings(nRows,1);
Statistic = strings(nRows,1);

row = 0;
for i = 1:nItems
    row = row + 1; Item(row) = ItemLabels(i); Statistic(row) = "mean";
    row = row + 1; Item(row) = "";           Statistic(row) = "std";
    row = row + 1; Item(row) = "";           Statistic(row) = "rank";
end
row = row + 1; Item(row) = "AvgRank"; Statistic(row) = "rank";

Tbl = table(Item, Statistic);
Tbl.Properties.VariableNames{1} = itemName;

for k = 1:nAlg
    a = idxKeep(k);
    alias = matlab.lang.makeValidName(char(Entries(a).alias));
    col = strings(nRows,1);
    row = 0;
    for i = 1:nItems
        row = row + 1; col(row) = sprintf('%.3e', MeanMat(i,a));
        row = row + 1; col(row) = sprintf('%.3e', StdMat(i,a));
        row = row + 1; col(row) = sprintf('%.3f', RankMat(i,a));
    end
    row = row + 1; col(row) = sprintf('%.3f', AvgRankAll(a));
    Tbl.(alias) = col;
end
end

function Tbl = build_anchor_stats_table_unified(Entries, Pairwise, WTL_global, AvgRankAll, anchorIdx)
n = numel(Entries);
aliases = string({Entries.alias})';
algnames = string({Entries.name})';
roles = string({Entries.role})';
pairNames = string(Pairwise.Algorithm);
wtlNames = string(WTL_global.Algorithm);

keep = false(n,1);
for i = 1:n
    if i == anchorIdx
        continue;
    end
    if roles(i) == "SARSA" || roles(i) == "RANDOM" || aliases(i) == "MPHB"
        keep(i) = true;
    end
end
idxKeep = find(keep);

Alias = strings(numel(idxKeep),1);
Algorithm = strings(numel(idxKeep),1);
Role = strings(numel(idxKeep),1);
GlobalWTL = strings(numel(idxKeep),1);
AvgRank = nan(numel(idxKeep),1);
RawP = nan(numel(idxKeep),1);
HolmP = nan(numel(idxKeep),1);
RejectHolm = false(numel(idxKeep),1);
BetterPct = nan(numel(idxKeep),1);

for j = 1:numel(idxKeep)
    i = idxKeep(j);
    Alias(j) = aliases(i);
    Algorithm(j) = algnames(i);
    Role(j) = roles(i);
    AvgRank(j) = AvgRankAll(i);

    pidx = find(pairNames == algnames(i), 1, 'first');
    if ~isempty(pidx)
        RawP(j) = Pairwise.RawP(pidx);
        HolmP(j) = Pairwise.HolmAdjP(pidx);
        RejectHolm(j) = Pairwise.RejectHolm(pidx);
        BetterPct(j) = Pairwise.BetterThanAnchorPct(pidx);
    end

    gidx = find(wtlNames == algnames(i), 1, 'first');
    if ~isempty(gidx)
        GlobalWTL(j) = string(WTL_global.WTL_Global(gidx));
    end
end

Tbl = table(Alias, Algorithm, Role, GlobalWTL, AvgRank, RawP, HolmP, RejectHolm, BetterPct);
Tbl = sortrows(Tbl, {'AvgRank','RawP'}, {'ascend','ascend'});
end

%% ========================================================================
% FIGURES
% ========================================================================
function export_convergence_panels_cec(RawCurve, RawEstFEs, Entries, CFG, fig_dir)
mkdir_if_needed(fig_dir);

labels = cell(numel(Entries),1);
for a = 1:numel(Entries)
    labels{a} = clean_plot_label(Entries(a).alias);
end

[~, ord] = sort(cellfun(@plot_sort_key, labels));
labelsOrd = labels(ord);

fig = figure('Visible','off','Color','w','Units','pixels','Position',[50 50 1520 900]);
t = tiledlayout(ceil(numel(CFG.Functions)/4), 4, 'TileSpacing','tight', 'Padding','compact');
legendHandles = gobjects(numel(ord),1);

for fi = 1:numel(CFG.Functions)
    ax = nexttile(t);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    for kk = 1:numel(ord)
        a = ord(kk);
        [xfe, ymean] = build_mean_curve_from_cells(RawCurve(a,fi,:), RawEstFEs(a,fi,:));
        if isempty(xfe), continue; end
        yplot = max(ymean(:).', eps);
        [mk, ls] = get_plot_style(kk);
        idxMarks = unique(round(linspace(1, numel(xfe), min(8, numel(xfe)))));
        h = plot(ax, xfe, yplot, 'LineWidth', 1.0, ...
            'LineStyle', ls, 'Marker', mk, 'MarkerIndices', idxMarks, ...
            'MarkerSize', 5);
        if fi == 1
            legendHandles(kk) = h;
        end
    end

    set(ax, 'XScale','linear', 'YScale','log');
    xlabel(ax, 'FEs', 'FontSize', 11);
    ylabel(ax, 'Mean error', 'FontSize', 11);
    title(ax, sprintf('F%d', CFG.Functions(fi)), 'Interpreter','none', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax, 'FontSize', 10);
end

valid = isgraphics(legendHandles);
lgd = legend(legendHandles(valid), labelsOrd(valid), ...
    'NumColumns', 8, 'Interpreter','none', 'Orientation', 'horizontal');
lgd.Box = 'off';
lgd.Layout.Tile = 'south';
lgd.FontSize = 13;
try
    lgd.ItemTokenSize = [28, 14];
catch
end
drawnow;
try
    lgd.Units = 'normalized';
    pos = lgd.Position;
    lgd.Position = [0.05, pos(2), 0.90, pos(4)];
catch
end

saveas(fig, fullfile(fig_dir, 'cec_convergence_panels_logFE.png'));
savefig(fig, fullfile(fig_dir, 'cec_convergence_panels_logFE.fig'));
close(fig);
end

%% ========================================================================
% LATEX WRITERS
% ========================================================================
function write_latex_alias_map_unified(T, fname)
fid = fopen(fname,'w'); assert(fid>0,'Could not create %s', fname);
fprintf(fid, '%% Auto-generated\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
fprintf(fid, '\\caption{Alias map used throughout the CEC2022 sensitivity and ablation analysis.}\n');
fprintf(fid, '\\begin{tabular}{l l l}\\hline\n');
fprintf(fid, 'Alias & Method & Definition \\\\ \\hline\n');
for i = 1:height(T)
    fprintf(fid, '%s & %s & %s \\\\ \n', ...
        escape_tex(T.Alias(i)), escape_tex(T.Method(i)), escape_tex(T.Definition(i)));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
end

function write_latex_sensitivity_panel_unified(Entries, idxKeep, ItemLabels, itemName, domainName, metricName, MeanMat, StdMat, RankMat, AvgRankAll, panelName, fname)
fid = fopen(fname,'w'); assert(fid>0,'Could not create %s', fname);

fprintf(fid, '%% Auto-generated\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
if strcmpi(panelName,'SARSA')
    cap = sprintf('%s sensitivity results for SARSA configurations (BEFORE + AFTER) with MPHB as reference. For each %s, the %s mean, standard deviation, and global rank are reported in separate rows. Global ranks are computed from the full method set using mean values. Lower values are better. The minimum mean within each displayed panel is highlighted in bold.', domainName, lower(itemName), metricName);
else
    cap = sprintf('%s sensitivity results for RANDOM configurations (BEFORE + AFTER) with MPHB as reference. For each %s, the %s mean, standard deviation, and global rank are reported in separate rows. Global ranks are computed from the full method set using mean values. Lower values are better. The minimum mean within each displayed panel is highlighted in bold.', domainName, lower(itemName), metricName);
end
fprintf(fid, '\\caption{%s}\n', escape_tex(cap));
fprintf(fid, '\\resizebox{\\textwidth}{!}{\n');
fprintf(fid, '\\begin{tabular}{ll');
for k = 1:numel(idxKeep)
    fprintf(fid, 'c');
end
fprintf(fid, '}\\n');
fprintf(fid, '\\hline\n');

fprintf(fid, '%s & Statistic', itemName);
for k = 1:numel(idxKeep)
    fprintf(fid, ' & %s', escape_tex(string(Entries(idxKeep(k)).alias)));
end
fprintf(fid, ' \\\\ \\hline\n');

for i = 1:numel(ItemLabels)
    vals = nan(1,numel(idxKeep));
    for k = 1:numel(idxKeep)
        vals(k) = MeanMat(i, idxKeep(k));
    end
    [~, bestk] = min(vals);

    fprintf(fid, '\\multirow{3}{*}{%s} & mean', escape_tex(ItemLabels(i)));
    for k = 1:numel(idxKeep)
        if k == bestk
            fprintf(fid, ' & \\textbf{%s}', sprintf('%.3e', MeanMat(i, idxKeep(k))));
        else
            fprintf(fid, ' & %s', sprintf('%.3e', MeanMat(i, idxKeep(k))));
        end
    end
    fprintf(fid, ' \\\\ \n');

    fprintf(fid, ' & std');
    for k = 1:numel(idxKeep)
        fprintf(fid, ' & %s', sprintf('%.3e', StdMat(i, idxKeep(k))));
    end
    fprintf(fid, ' \\\\ \n');

    fprintf(fid, ' & rank');
    for k = 1:numel(idxKeep)
        fprintf(fid, ' & %s', sprintf('%.3f', RankMat(i, idxKeep(k))));
    end
    fprintf(fid, ' \\\\ \\hline\n');
end

avg = nan(1,numel(idxKeep));
for k = 1:numel(idxKeep), avg(k) = AvgRankAll(idxKeep(k)); end
[~, bestk] = min(avg);
fprintf(fid, 'AvgRank & rank');
for k = 1:numel(idxKeep)
    if k == bestk
        fprintf(fid, ' & \\textbf{%s}', sprintf('%.3f', AvgRankAll(idxKeep(k))));
    else
        fprintf(fid, ' & %s', sprintf('%.3f', AvgRankAll(idxKeep(k))));
    end
end
fprintf(fid, ' \\\\ \n');

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
end

function write_latex_anchor_stats_unified(T, anchorAlias, fname)
fid = fopen(fname,'w'); assert(fid>0,'Could not create %s', fname);
fprintf(fid, '%% Auto-generated\n');
fprintf(fid, '\\begin{table*}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
cap = sprintf('Global statistical comparison against the best SARSA configuration (anchor = %s). W/T/L is reported from the perspective of the competitor. Lower AvgRank is better.', char(anchorAlias));
fprintf(fid, '\\caption{%s}\n', escape_tex(cap));
fprintf(fid, '\\begin{tabular}{l l c c c c c c}\\hline\n');
fprintf(fid, 'Method & Role & Global W/T/L & AvgRank & $p$ & Holm $p$ & Reject & Better(\\%%) \\\\ \\hline\n');
for i = 1:height(T)
    fprintf(fid, '%s & %s & %s & %.3f & %s & %s & %d & %.2f \\\\ \n', ...
        escape_tex(T.Alias(i)), escape_tex(T.Role(i)), escape_tex(T.GlobalWTL(i)), T.AvgRank(i), ...
        fmt_p(T.RawP(i)), fmt_p(T.HolmP(i)), T.RejectHolm(i), T.BetterPct(i));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table*}\n');
fclose(fid);
end

%% ========================================================================
% GENERIC HELPERS
% ========================================================================
function [xFE, yMean] = build_mean_curve_from_cells(curveCellRow, feRow)
n = numel(curveCellRow);
curves = {};
fes = [];
for i = 1:n
    c = curveCellRow{i};
    if isempty(c) || ~isnumeric(c), continue; end
    c = c(:).';
    if all(isnan(c)), continue; end
    curves{end+1,1} = c; %#ok<AGROW>
    fe = feRow(1,1,i);
    if isfinite(fe) && fe > 0, fes(end+1,1) = fe; end %#ok<AGROW>
end
if isempty(curves)
    xFE = [];
    yMean = [];
    return;
end
L = min(cellfun(@numel, curves));
M = nan(numel(curves), L);
for i = 1:numel(curves)
    M(i,:) = curves{i}(1:L);
end
yMean = mean(M, 1, 'omitnan');
feEnd = mean(fes, 'omitnan');
if ~isfinite(feEnd) || feEnd <= 0, feEnd = L; end
xFE = linspace(max(feEnd/L,1), feEnd, L);
end

function [mk, ls] = get_plot_style(k)
markers = {'o','s','d','^','v','>','<','p','h','x','+','*','.','|','_'};
linestyles = {'-','--',':','-.'};
nM = numel(markers);
nL = numel(linestyles);
idx = k - 1;
mk = markers{mod(idx, nM) + 1};
ls = linestyles{mod(floor(idx / nM), nL) + 1};
end

function lbl = clean_plot_label(alias)
alias = string(alias);
aliasCompact = replace(alias, "-", "");
if alias == "MPHB"
    lbl = 'MPHB';
elseif startsWith(alias, "RB-") || startsWith(alias, "RA-")
    lbl = "MPHBR-" + aliasCompact;
elseif startsWith(alias, "B-") || startsWith(alias, "A-")
    lbl = "MPHBS-" + aliasCompact;
else
    lbl = aliasCompact;
end
lbl = char(lbl);
end

function key = plot_sort_key(lbl)
if strcmp(lbl,'MPHB')
    key = 1e9; return;
end
if startsWith(lbl,'MPHBS')
    base = 0;
elseif startsWith(lbl,'MPHBR')
    base = 1e5;
else
    base = 2e5;
end
nums = regexp(lbl, 'C(\d+)', 'tokens', 'once');
cfg = 99;
if ~isempty(nums), cfg = str2double(nums{1}); end
if contains(lbl, '-B-')
    fad = 0;
elseif contains(lbl, '-A-')
    fad = 1;
else
    fad = 2;
end
key = base + fad * 100 + cfg;
end

function s = fmt_p(x)
if isnan(x), s = 'NaN'; else, s = sprintf('%.3e', x); end
end

function s = escape_tex(x)
if isstring(x) || ischar(x)
    s = char(x);
else
    s = char(string(x));
end
s = strrep(s, '\', '\textbackslash{}');
s = strrep(s, '_', '\_');
s = strrep(s, '%', '\%');
s = strrep(s, '&', '\&');
s = strrep(s, '#', '\#');
s = strrep(s, '{', '\{');
s = strrep(s, '}', '\}');
s = strrep(s, '^', '\^{}');
s = strrep(s, '~', '\~{}');
end

function mkdir_if_needed(p)
if exist(p, 'dir') ~= 7
    mkdir(p);
end
end

function r = tiedrank_local(x)
x = x(:);
r = nan(size(x));
mask = ~isnan(x);
if ~any(mask), return; end
z = x(mask);
[sortedX, ord] = sort(z);
rk = nan(size(z));
i = 1;
while i <= numel(z)
    j = i;
    while j < numel(z) && isequaln(sortedX(j+1), sortedX(i))
        j = j + 1;
    end
    rankVal = (i + j) / 2;
    rk(ord(i:j)) = rankVal;
    i = j + 1;
end
r(mask) = rk;
end

function [p, h, stats] = signrank_safe(x, y)
mask = ~(isnan(x) | isnan(y));
x = x(mask); y = y(mask);
if isempty(x)
    p = nan; h = false; stats = struct(); return;
end
try
    [p, h, stats] = signrank(x, y, 'method', 'approximate');
catch
    try
        [p, h, stats] = signrank(x, y);
    catch
        p = nan; h = false; stats = struct();
    end
end
end

function [adjp, reject] = holm_correction(pvals, alpha)
pvals = pvals(:);
m = numel(pvals);
[ps, ord] = sort(pvals, 'ascend');

adj_sorted = nan(m,1);
for i = 1:m
    adj_sorted(i) = (m - i + 1) * ps(i);
end
adj_sorted = min(adj_sorted, 1);
for i = 2:m
    adj_sorted(i) = max(adj_sorted(i), adj_sorted(i-1));
end

reject_sorted = false(m,1);
for i = 1:m
    thresh = alpha / (m - i + 1);
    if ps(i) <= thresh
        reject_sorted(i) = true;
    else
        reject_sorted(i:end) = false;
        break;
    end
end

adjp = nan(m,1);
reject = false(m,1);
adjp(ord) = adj_sorted;
reject(ord) = reject_sorted;
end
