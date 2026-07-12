function run_FIR_ablation_AC1_common_seeds
% =========================================================================
% RUN_FIR_SENSITIVITY_COMMON_SEEDS
% =========================================================================
% FIR ablation for the best-ranked SARSA sensitivity setting (A-C1).
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
%     fir_alias_map
%     fir_sarsa_sensitivity
%     fir_random_sensitivity
%     fir_anchor_stats
%
%   Figures:
%     fir_convergence_panels_logFE.png / .fig
%
% Design principles
% -----------------
% - Rank for each FIR case is computed from the mean fitness of the FULL
%   method set.
% - AvgRank is the mean of the per-case global ranks.
% - Panel tables (SARSA / RANDOM) only FILTER displayed columns; they do
%   NOT recompute ranks on the displayed subset.
% - Anchor is selected among SARSA methods using global AvgRank, breaking
%   ties by GlobalMeanFit, to mirror the CEC reference logic.
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
CFG.EXP_NAME      = ['FIR_ablation_AC1_' CFG.VERSION_TAG];
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
mkdir_if_needed(fullfile(CFG.OUTDIR, 'checkpoints_tasks'));
mkdir_if_needed(fullfile(CFG.OUTDIR, 'checkpoints_cases'));

CFG.SELF_COPY = true;
if CFG.SELF_COPY
    try
        thisFile = mfilename('fullpath');
        copyfile([thisFile '.m'], fullfile(CFG.OUTDIR, 'logs', [mfilename '.m']));
    catch
    end
end

%% ========================= GLOBAL SETTINGS =========================
CFG.PopSize = 30;
CFG.nRuns   = 30;
CFG.Max_FEs = 150000;
CFG.LB      = -1;
CFG.UB      = 1;
CFG.N_fft   = 2048;
CFG.case_list = 1:8;

CFG.CALIBRATE_FE_MATCH = true;
CFG.T_PROBE1 = 1;
CFG.T_PROBE2 = 2;

CFG.USE_PARFOR = true;
CFG.N_WORKERS  = 8;
CFG.EXPORT_CSV = true;
CFG.EXPORT_TEX = true;
CFG.EXPORT_FIG = true;
CFG.SAVE_RAW_MAT = true;

CFG.MASTER_SEED = 20260405;
CFG.TIE_TOL     = 1e-12;
CFG.RESUME_IF_AVAILABLE = true;
CFG.RETRY_FAILED_TASKS  = true;
CFG.PRINT_CASE_PROGRESS = true;
CFG.CLOSE_POOL_AT_END   = true;

CFG.FIXED_K  = 3;
CFG.FIXED_NI = 9;

%% ========================= ABLATION SELECTION =========================
% Best-ranked SARSA setting and its matched RANDOM mirror.
CFG.ANCHOR_CFG_ID          = 'C1';
CFG.ANCHOR_FAD             = 'AFTER';
CFG.MATCHED_RANDOM_CFG_ID  = 'C1';
CFG.MATCHED_RANDOM_FAD     = 'AFTER';
CFG.BEST_RANDOM_CFG_ID     = 'C2';
CFG.BEST_RANDOM_FAD        = 'BEFORE';

USER_CONFIGS = struct([]);

%% ========================= FIR CASES =========================
Cases = build_fir_cases(CFG.N_fft);
Cases = Cases(CFG.case_list);
numCases = numel(Cases);
CFG.PROBE_DIM = Cases(1).Dim;

%% ========================= SANITY CHECKS =========================
needFiles = {'MPHBS_main', ...
             'MPHBS_random_mirror', ...
             'MPHB_baseline', 'HBA', 'MPA'};
for k = 1:numel(needFiles)
    if isempty(which(needFiles{k}))
        error('Required function not found on MATLAB path: %s', needFiles{k});
    end
end

%% ========================= ENTRIES =========================
Entries = build_entries_full(USER_CONFIGS, CFG);
numAlgs = numel(Entries);
for a = 1:numAlgs
    [Entries(a).T, Entries(a).EstimatedFEs] = estimate_T_and_FEs(Entries(a), CFG);
end

save(fullfile(CFG.OUTDIR, 'mat', 'config_snapshot.mat'), 'CFG', 'USER_CONFIGS', 'Entries', 'Cases');
write_entries_manifest(Entries, fullfile(CFG.OUTDIR, 'csv', 'entries_manifest.csv'));

fprintf('\n============================================================\n');
fprintf('RUNNING %s\n', CFG.EXP_NAME);
fprintf('Algorithms/entries: %d\n', numAlgs);
fprintf('Cases: %d | Runs per case: %d\n', numCases, CFG.nRuns);
fprintf('PopSize=%d | Max_FEs=%d | Workers=%d\n', CFG.PopSize, CFG.Max_FEs, CFG.N_WORKERS);
fprintf('Fixed K=%d | Fixed Ni=%d\n', CFG.FIXED_K, CFG.FIXED_NI);
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
RawFinal   = nan(numAlgs, numCases, CFG.nRuns);
TimeAlg    = nan(numAlgs, numCases, CFG.nRuns);
FEAlg      = nan(numAlgs, numCases, CFG.nRuns);
CurvesRaw  = cell(numAlgs, numCases, CFG.nRuns);
BestPos    = cell(numAlgs, numCases, CFG.nRuns);
DoneMask   = false(numAlgs, numCases, CFG.nRuns);
OKMask     = false(numAlgs, numCases, CFG.nRuns);
MsgMask    = strings(numAlgs, numCases, CFG.nRuns);
SeedMask   = nan(numAlgs, numCases, CFG.nRuns);

%% ========================= RESUME CASE CHECKPOINTS =========================
for ci = 1:numCases
    Case = Cases(ci);
    case_ckpt_file = fullfile(CFG.OUTDIR, 'checkpoints_cases', sprintf('C%02d_done.mat', Case.id));
    if exist(case_ckpt_file, 'file') == 2 && CFG.RESUME_IF_AVAILABLE
        S = load(case_ckpt_file, 'caseData');
        if isfield(S, 'caseData')
            [RawFinal, TimeAlg, FEAlg, CurvesRaw, BestPos, DoneMask, OKMask, MsgMask, SeedMask] = ...
                merge_case_checkpoint(S.caseData, ci, RawFinal, TimeAlg, FEAlg, CurvesRaw, BestPos, DoneMask, OKMask, MsgMask, SeedMask);
        end
    end
end

%% ========================= EVALUATION BY CASE =========================
for ci = 1:numCases
    Case = Cases(ci);
    fprintf('------------------------------------------------------------\n');
    fprintf('Case C%02d | %s | %s | Dim=%d\n', Case.id, Case.name, Case.type, Case.Dim);

    case_ckpt_file = fullfile(CFG.OUTDIR, 'checkpoints_cases', sprintf('C%02d_done.mat', Case.id));
    done_ci = squeeze(DoneMask(:,ci,:));
    missing = find(~done_ci);

    if CFG.PRINT_CASE_PROGRESS
        fprintf('Completed tasks found in checkpoint: %d / %d\n', nnz(done_ci), numAlgs * CFG.nRuns);
    end

    if isempty(missing)
        fprintf('C%02d already complete. Loading checkpoint only.\n', Case.id);
        continue;
    end

    LocalTaskTable = build_local_task_table(numAlgs, CFG.nRuns, missing);
    numLocalTasks = size(LocalTaskTable, 1);
    ResultsCell = cell(numLocalTasks, 1);

    if CFG.PRINT_CASE_PROGRESS
        fprintf('C%02d missing tasks to execute: %d\n', Case.id, numLocalTasks);
    end

    if CFG.USE_PARFOR
        parfor taskId = 1:numLocalTasks
            taskFile = build_task_checkpoint_file(fullfile(CFG.OUTDIR,'checkpoints_tasks'), Case.id, LocalTaskTable(taskId,2), Entries(LocalTaskTable(taskId,1)).name);
            ResultsCell{taskId} = load_or_run_single_task_fir(taskFile, Case, LocalTaskTable(taskId,2), Entries(LocalTaskTable(taskId,1)), CFG);
        end
    else
        for taskId = 1:numLocalTasks
            taskFile = build_task_checkpoint_file(fullfile(CFG.OUTDIR,'checkpoints_tasks'), Case.id, LocalTaskTable(taskId,2), Entries(LocalTaskTable(taskId,1)).name);
            ResultsCell{taskId} = load_or_run_single_task_fir(taskFile, Case, LocalTaskTable(taskId,2), Entries(LocalTaskTable(taskId,1)), CFG);
        end
    end

    raw_ci    = squeeze(RawFinal(:,ci,:));
    time_ci   = squeeze(TimeAlg(:,ci,:));
    fe_ci     = squeeze(FEAlg(:,ci,:));
    curves_ci = squeeze(CurvesRaw(:,ci,:));
    best_ci   = squeeze(BestPos(:,ci,:));
    done_ci   = squeeze(DoneMask(:,ci,:));
    ok_ci     = squeeze(OKMask(:,ci,:));
    msg_ci    = squeeze(MsgMask(:,ci,:));
    seed_ci   = squeeze(SeedMask(:,ci,:));

    for taskId = 1:numLocalTasks
        rr = ResultsCell{taskId};
        a = rr.algIdx;
        r = rr.runIdx;
        raw_ci(a,r)    = rr.best_f;
        time_ci(a,r)   = rr.elapsed;
        fe_ci(a,r)     = rr.fe_used;
        curves_ci{a,r} = rr.curve;
        best_ci{a,r}   = rr.best_x;
        done_ci(a,r)   = rr.ok;
        ok_ci(a,r)     = rr.ok;
        msg_ci(a,r)    = string(rr.error_message);
        seed_ci(a,r)   = rr.seed;
    end

    RawFinal(:,ci,:) = raw_ci;
    TimeAlg(:,ci,:)  = time_ci;
    FEAlg(:,ci,:)    = fe_ci;
    CurvesRaw(:,ci,:)= curves_ci;
    BestPos(:,ci,:)  = best_ci;
    DoneMask(:,ci,:) = done_ci;
    OKMask(:,ci,:)   = ok_ci;
    MsgMask(:,ci,:)  = msg_ci;
    SeedMask(:,ci,:) = seed_ci;

    caseData = struct();
    caseData.raw_ci    = raw_ci;
    caseData.time_ci   = time_ci;
    caseData.fe_ci     = fe_ci;
    caseData.curves_ci = curves_ci;
    caseData.best_ci   = best_ci;
    caseData.done_ci   = done_ci;
    caseData.ok_ci     = ok_ci;
    caseData.msg_ci    = msg_ci;
    caseData.seed_ci   = seed_ci;
    save(case_ckpt_file, 'caseData', '-v7.3');

    if CFG.PRINT_CASE_PROGRESS
        fprintf('C%02d checkpoint saved: %d / %d tasks complete\n', Case.id, nnz(done_ci), numAlgs * CFG.nRuns);
    end
end

%% ========================= TASK RESULTS =========================
TaskResults = rebuild_task_results_from_arrays(RawFinal, TimeAlg, FEAlg, SeedMask, OKMask, MsgMask, Entries, Cases);
if CFG.EXPORT_CSV
    writetable(TaskResults, fullfile(CFG.OUTDIR, 'csv', 'task_results.csv'));
end

%% ========================= SUMMARY / STATS =========================
RanksRunBlock = nan(numAlgs, numCases, CFG.nRuns);
for ci = 1:numCases
    for r = 1:CFG.nRuns
        vec = RawFinal(:,ci,r);
        if all(isnan(vec)), continue; end
        RanksRunBlock(:,ci,r) = tiedrank_local(vec')';
    end
end

Summary = build_summary_table(Entries, RawFinal, TimeAlg, FEAlg, OKMask, DoneMask, RanksRunBlock);
[Y_valid, p_friedman, MeanRanksFromFriedman] = build_friedman_from_fits(RawFinal);
for a = 1:numAlgs
    Summary.FriedmanMeanRank(a) = MeanRanksFromFriedman(a);
end

[MeanByCaseAll, StdByCaseAll, RankByMeanCaseAll, AvgRankAll] = build_casewise_stats_from_mean(RawFinal);
Summary.AvgRankAll = AvgRankAll;
Summary = sortrows(Summary, 'FriedmanMeanRank', 'ascend');

anchorAliasTarget = build_alias('ANCHOR', CFG.ANCHOR_FAD, CFG.ANCHOR_CFG_ID);
anchorIdx = find(strcmpi(string({Entries.alias}), string(anchorAliasTarget)), 1, 'first');
if isempty(anchorIdx), error('Could not find requested anchor alias: %s', anchorAliasTarget); end
anchorName = Entries(anchorIdx).name;
anchorAlias = Entries(anchorIdx).alias;

Pairwise      = build_pairwise_vs_anchor(Entries, RawFinal, anchorIdx);
WTL_by_case   = build_wtl_by_case(RawFinal, Entries, anchorIdx, Cases, CFG.TIE_TOL);
WTL_global    = build_wtl_global(RawFinal, Entries, anchorIdx, CFG.TIE_TOL);

%% ========================= BASE EXPORTS =========================
if CFG.EXPORT_CSV
    writetable(Summary,       fullfile(CFG.OUTDIR, 'csv', 'summary_global.csv'));
    writetable(Pairwise,      fullfile(CFG.OUTDIR, 'csv', 'pairwise_vs_anchor.csv'));
    writetable(WTL_by_case,   fullfile(CFG.OUTDIR, 'csv', 'wtl_by_case_vs_anchor.csv'));
    writetable(WTL_global,    fullfile(CFG.OUTDIR, 'csv', 'wtl_global_vs_anchor.csv'));
end

%% ========================= MAIN-TEXT TABLES =========================
export_fir_ablation_tables(Entries, Pairwise, WTL_global, Cases, ...
    anchorIdx, string(anchorAlias), MeanByCaseAll, StdByCaseAll, RankByMeanCaseAll, AvgRankAll, ...
    fullfile(CFG.OUTDIR,'tables'), fullfile(CFG.OUTDIR,'csv'));

%% ========================= FIGURES =========================
if CFG.EXPORT_FIG
    export_convergence_panels_fir(CurvesRaw, FEAlg, Entries, Cases, fullfile(CFG.OUTDIR, 'figures'));
end

%% ========================= FINAL SAVE =========================
if CFG.SAVE_RAW_MAT
    save(fullfile(CFG.OUTDIR, 'mat', 'final_workspace.mat'), ...
        'CFG', 'USER_CONFIGS', 'Entries', 'Cases', 'Summary', 'Pairwise', 'WTL_by_case', 'WTL_global', ...
        'TaskResults', 'RawFinal', 'TimeAlg', 'FEAlg', 'CurvesRaw', 'BestPos', 'DoneMask', 'OKMask', 'MsgMask', 'SeedMask', ...
        'RanksRunBlock', 'MeanByCaseAll', 'StdByCaseAll', 'RankByMeanCaseAll', 'AvgRankAll', ...
        'Y_valid', 'p_friedman', 'anchorIdx', 'anchorName', 'anchorAlias', '-v7.3');
end

fprintf('\n============================================================\n');
fprintf('%s DONE\n', CFG.EXP_NAME);
fprintf('Anchor for pairwise tests: %s (%s)\n', anchorAlias, anchorName);
fprintf('Friedman p-value = %.6g\n', p_friedman);
disp(Summary(:, {'Alias','Algorithm','Family','Role','GlobalMeanFit','GlobalStdFit','GlobalMedianFit','GlobalMeanTimeSec','GlobalMeanFEs','AvgRankAll','FriedmanMeanRank'}));
fprintf('Outputs saved in: %s\n', CFG.OUTDIR);
fprintf('============================================================\n');

if CFG.CLOSE_POOL_AT_END
    pool = gcp('nocreate');
    if ~isempty(pool), delete(pool); end
end
end

%% ========================================================================
% FIR CASES / OBJECTIVE
% ========================================================================
function Cases = build_fir_cases(N_fft)
Cases(1) = make_case(1, 'LPF_Easy',    'LPF', 30, N_fft, 0.18, 0.25, [],   [],   1, 100);
Cases(2) = make_case(2, 'LPF_Hard',    'LPF', 30, N_fft, 0.22, 0.26, [],   [],   1, 120);
Cases(3) = make_case(3, 'HPF_Easy',    'HPF', 30, N_fft, [],   [],   0.75, 0.82, 1, 100);
Cases(4) = make_case(4, 'HPF_Hard',    'HPF', 30, N_fft, [],   [],   0.74, 0.78, 1, 120);
Cases(5) = make_case(5, 'BPF_Wide',    'BPF', 30, N_fft, 0.25, 0.35, 0.65, 0.75, 1, 100);
Cases(6) = make_case(6, 'BPF_Narrow',  'BPF', 30, N_fft, 0.32, 0.40, 0.60, 0.68, 1, 120);
Cases(7) = make_case(7, 'BSF_Wide',    'BSF', 30, N_fft, 0.25, 0.35, 0.65, 0.75, 1, 100);
Cases(8) = make_case(8, 'BSF_Narrow',  'BSF', 30, N_fft, 0.32, 0.40, 0.60, 0.68, 1, 120);
end

function C = make_case(id, name, type, N_order, N_fft, f1, f2, f3, f4, Wp, Ws)
Dim = N_order + 1;
w_norm = (0:(N_fft/2)) / (N_fft/2);
switch upper(type)
    case 'LPF'
        Fp = f1; Fs = f2;
        passband_idx = find(w_norm <= Fp);
        stopband_idx = find(w_norm >= Fs);
    case 'HPF'
        Fs = f3; Fp = f4;
        passband_idx = find(w_norm >= Fp);
        stopband_idx = find(w_norm <= Fs);
    case 'BPF'
        Fs_stop1 = f1; Fp_pass1 = f2; Fp_pass2 = f3; Fs_stop2 = f4;
        passband_idx = find(w_norm >= Fp_pass1 & w_norm <= Fp_pass2);
        stopband_idx = find(w_norm <= Fs_stop1 | w_norm >= Fs_stop2);
    case 'BSF'
        Fp_pass1 = f1; Fs_stop1 = f2; Fs_stop2 = f3; Fp_pass2 = f4;
        passband_idx = find(w_norm <= Fp_pass1 | w_norm >= Fp_pass2);
        stopband_idx = find(w_norm >= Fs_stop1 & w_norm <= Fs_stop2);
    otherwise
        error('Unsupported FIR type: %s', type);
end
C = struct();
C.id = id; C.name = name; C.type = upper(type); C.N_order = N_order; C.Dim = N_order + 1; C.N_fft = N_fft;
C.w_norm = w_norm; C.passband_idx = passband_idx; C.stopband_idx = stopband_idx;
C.inv_Lp = 1 / max(1, numel(passband_idx));
C.inv_Ls = 1 / max(1, numel(stopband_idx));
C.Wp = Wp; C.Ws = Ws;
C.spec = struct('f1',f1,'f2',f2,'f3',f3,'f4',f4);
end

function y = smartObjectiveFIR(x, f_raw)
global FE_COUNTER_FIR
[r, c] = size(x);
if r > 1 && c > 1
    y = zeros(c, 1);
    for i = 1:c
        FE_COUNTER_FIR = FE_COUNTER_FIR + 1;
        indiv = x(:, i);
        y(i) = f_raw(indiv(:).');
    end
else
    FE_COUNTER_FIR = FE_COUNTER_FIR + 1;
    y = f_raw(x(:).');
end
end

function fitness = fir_fitness_generic(coeffs, N_fft, idx_stop, idx_pass, inv_Lp, inv_Ls, Wp, Ws)
H_fft = fft(coeffs, N_fft);
H_half = H_fft(1:(N_fft/2 + 1));
if ~isempty(idx_stop)
    H_stop = H_half(idx_stop);
    MSEs = inv_Ls * sum(real(H_stop).^2 + imag(H_stop).^2);
else
    MSEs = 0;
end
if ~isempty(idx_pass)
    H_pass_mag = abs(H_half(idx_pass));
    MSEp = inv_Lp * sum((H_pass_mag - 1).^2);
else
    MSEp = 0;
end
fitness = Wp * MSEp + Ws * MSEs;
end

%% ========================================================================
% TASK EXECUTION
% ========================================================================
function res = load_or_run_single_task_fir(taskFile, Case, runIdx, Entry, CFG)
if exist(taskFile, 'file') == 2
    S = load(taskFile, 'res');
    if isfield(S,'res') && isfield(S.res,'ok') && S.res.ok
        res = S.res;
        return;
    end
end

% Common random numbers across ablation methods for this case/run block.
seed_now = CFG.MASTER_SEED + 1000 * Case.id + runIdx;
rng(seed_now, 'twister');

global FE_COUNTER_FIR
FE_COUNTER_FIR = 0;
f_raw = @(coeffs) fir_fitness_generic(coeffs, Case.N_fft, Case.stopband_idx, Case.passband_idx, Case.inv_Lp, Case.inv_Ls, Case.Wp, Case.Ws);
F = @(x) smartObjectiveFIR(x, f_raw);

t0 = tic;
res = struct();
res.algIdx = Entry.alg_index;
res.runIdx = runIdx;
res.case_id = Case.id;
res.seed = seed_now;
res.ok = false;
res.error_message = "";

try
    switch upper(Entry.family)
        case 'MPHBS_SARSA'
            [bestVal, bestX, conv] = MPHBS_main(CFG.PopSize, Entry.T, CFG.LB, CFG.UB, Case.Dim, F, ...
                Entry.K, Entry.Ni, Entry.rho_sub, Entry.w_best, Entry.gamma, Entry.fad_before_p2, Entry.fad_after_p2);
        case 'MPHBS_RANDOM'
            [bestVal, bestX, conv] = MPHBS_random_mirror(CFG.PopSize, Entry.T, CFG.LB, CFG.UB, Case.Dim, F, ...
                Entry.K, Entry.Ni, Entry.rho_sub, Entry.w_best, Entry.gamma, Entry.fad_before_p2, Entry.fad_after_p2);
        case 'MPHB'
            [bestVal, bestX, conv] = MPHB_baseline(CFG.PopSize, Entry.T, CFG.LB, CFG.UB, Case.Dim, F);
        case 'HBA_BASE'
            [bestX, bestVal, conv] = HBA(F, Case.Dim, CFG.LB, CFG.UB, Entry.T, CFG.PopSize);
        case 'MPA_BASE'
            [bestVal, bestX, conv] = MPA(CFG.PopSize, Entry.T, CFG.LB, CFG.UB, Case.Dim, F);
        otherwise
            error('Unknown family: %s', Entry.family);
    end
    elapsed = toc(t0);
    res.ok = true;
    res.best_f = bestVal;
    res.best_x = bestX;
    res.elapsed = elapsed;
    res.fe_used = FE_COUNTER_FIR;
    res.curve = conv(:).';
catch ME
    elapsed = toc(t0);
    res.ok = false;
    res.best_f = NaN;
    res.best_x = [];
    res.elapsed = elapsed;
    res.fe_used = FE_COUNTER_FIR;
    res.curve = [];
    res.error_message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
end
save(taskFile, 'res', '-v7.3');
end

%% ========================================================================
% ENTRY BUILDERS
% ========================================================================
function cfg = make_cfg(cfg_id, w_best, gamma, rho_sub, fad_name, enabled)
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

function Entries = build_entries_full(USER_CONFIGS, CFG) %#ok<INUSD>
Entries = struct('pair_id', {}, 'cfg_id', {}, 'name', {}, 'family', {}, 'role', {}, 'K', {}, 'Ni', {}, ...
    'gamma', {}, 'w_best', {}, 'rho_sub', {}, 'label', {}, 'fad_name', {}, 'fad_before_p2', {}, 'fad_after_p2', {}, ...
    'alg_index', {}, 'T', {}, 'EstimatedFEs', {}, 'alias', {});
idx = 0;
idx = idx + 1;
Entries(idx) = make_entry(CFG.ANCHOR_CFG_ID, "ANCHOR_SARSA", 'MPHBS_SARSA', 'ANCHOR', 'anchor_sarsa', ...
    ablation_cfg_param('K', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.ANCHOR_CFG_ID, CFG.ANCHOR_FAD, CFG), ...
    CFG.ANCHOR_FAD, strcmpi(CFG.ANCHOR_FAD,'BEFORE'), strcmpi(CFG.ANCHOR_FAD,'AFTER'), idx, build_alias('ANCHOR', CFG.ANCHOR_FAD, CFG.ANCHOR_CFG_ID));
idx = idx + 1;
Entries(idx) = make_entry(CFG.MATCHED_RANDOM_CFG_ID, "MATCHED_RANDOM", 'MPHBS_RANDOM', 'MATCHED_RANDOM', 'matched_random', ...
    ablation_cfg_param('K', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.MATCHED_RANDOM_CFG_ID, CFG.MATCHED_RANDOM_FAD, CFG), ...
    CFG.MATCHED_RANDOM_FAD, strcmpi(CFG.MATCHED_RANDOM_FAD,'BEFORE'), strcmpi(CFG.MATCHED_RANDOM_FAD,'AFTER'), idx, build_alias('MATCHED_RANDOM', CFG.MATCHED_RANDOM_FAD, CFG.MATCHED_RANDOM_CFG_ID));
idx = idx + 1;
Entries(idx) = make_entry(CFG.BEST_RANDOM_CFG_ID, "BEST_RANDOM", 'MPHBS_RANDOM', 'BEST_RANDOM', 'best_random', ...
    ablation_cfg_param('K', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('Ni', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('gamma', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('w_best', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    ablation_cfg_param('rho_sub', CFG.BEST_RANDOM_CFG_ID, CFG.BEST_RANDOM_FAD, CFG), ...
    CFG.BEST_RANDOM_FAD, strcmpi(CFG.BEST_RANDOM_FAD,'BEFORE'), strcmpi(CFG.BEST_RANDOM_FAD,'AFTER'), idx, 'R-best');
idx = idx + 1;
Entries(idx) = make_entry('BASE', 'BASE', 'MPHB', 'BASELINE', 'baseline', nan, nan, nan, nan, nan, 'BASE', false, false, idx, 'MPHB');
idx = idx + 1;
Entries(idx) = make_entry('HBA', 'HBA', 'HBA_BASE', 'PARENT', 'parent_hba', nan, nan, nan, nan, nan, 'BASE', false, false, idx, 'HBA');
idx = idx + 1;
Entries(idx) = make_entry('MPA', 'MPA', 'MPA_BASE', 'PARENT', 'parent_mpa', nan, nan, nan, nan, nan, 'BASE', false, false, idx, 'MPA');
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


function e = make_entry(cfg_id, pair_id, family, role, label, K, Ni, gamma, w_best, rho_sub, fad_name, fad_before_p2, fad_after_p2, alg_index, alias)
e = struct();
e.pair_id = char(pair_id);
e.cfg_id  = char(cfg_id);
e.family  = char(family);
e.role    = char(role);
e.K       = K;
e.Ni      = Ni;
e.gamma   = gamma;
e.w_best  = w_best;
e.rho_sub = rho_sub;
e.fad_name = char(fad_name);
e.fad_before_p2 = logical(fad_before_p2);
e.fad_after_p2  = logical(fad_after_p2);
e.label   = char(label);
e.alg_index = alg_index;
e.T = nan;
e.EstimatedFEs = nan;
e.alias = char(alias);
e.name = build_entry_name(e);
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
        name = sprintf('MPHBS_SARSA_%s_K%02d_Ni%02d_%s_W%03d_G%03d_R%03d', e.cfg_id, e.K, e.Ni, upper(e.fad_name), round(1000*e.w_best), round(1000*e.gamma), round(1000*e.rho_sub));
    case 'MPHBS_RANDOM'
        name = sprintf('MPHBS_RANDOM_%s_K%02d_Ni%02d_%s_W%03d_G%03d_R%03d', e.cfg_id, e.K, e.Ni, upper(e.fad_name), round(1000*e.w_best), round(1000*e.gamma), round(1000*e.rho_sub));
    case 'MPHB'
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
MaxFE = CFG.Max_FEs;
probeDim = CFG.PROBE_DIM;
T1 = CFG.T_PROBE1;
T2 = CFG.T_PROBE2;
fe1 = probe_actual_fes(entry, CFG, probeDim, T1);
fe2 = probe_actual_fes(entry, CFG, probeDim, T2);
slope = fe2 - fe1;
intercept = fe1 - slope * T1;
if ~(isfinite(slope) && isfinite(intercept)) || slope <= 0
    error('Invalid FE probe for family=%s, cfg=%s (fe1=%g, fe2=%g).', entry.family, entry.cfg_id, fe1, fe2);
end
T = floor((MaxFE - intercept) / slope);
T = max(1, T);
estFE = intercept + slope * T;
while estFE > MaxFE && T > 1
    T = T - 1;
    estFE = intercept + slope * T;
end
end

function fe_used = probe_actual_fes(entry, CFG, probeDim, Tprobe)
global FE_COUNTER_FIR
FE_COUNTER_FIR = 0;
rng(12345, 'twister');
f_raw = @(x) sum(x(:).^2);
F = @(x) probe_count_objective(x, f_raw);
switch upper(entry.family)
    case 'MPHBS_SARSA'
        MPHBS_main(CFG.PopSize, Tprobe, CFG.LB, CFG.UB, probeDim, F, entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, entry.fad_before_p2, entry.fad_after_p2);
    case 'MPHBS_RANDOM'
        MPHBS_random_mirror(CFG.PopSize, Tprobe, CFG.LB, CFG.UB, probeDim, F, entry.K, entry.Ni, entry.rho_sub, entry.w_best, entry.gamma, entry.fad_before_p2, entry.fad_after_p2);
    case 'MPHB'
        MPHB_baseline(CFG.PopSize, Tprobe, CFG.LB, CFG.UB, probeDim, F);
    case 'HBA_BASE'
        HBA(F, probeDim, CFG.LB, CFG.UB, Tprobe, CFG.PopSize);
    case 'MPA_BASE'
        MPA(CFG.PopSize, Tprobe, CFG.LB, CFG.UB, probeDim, F);
    otherwise
        error('Unknown family for FE probe: %s', entry.family);
end
fe_used = FE_COUNTER_FIR;
end

function y = probe_count_objective(x, f_raw)
global FE_COUNTER_FIR
[r, c] = size(x);
if r > 1 && c > 1
    y = zeros(c, 1);
    for i = 1:c
        FE_COUNTER_FIR = FE_COUNTER_FIR + 1;
        indiv = x(:, i);
        y(i) = f_raw(indiv(:).');
    end
else
    FE_COUNTER_FIR = FE_COUNTER_FIR + 1;
    y = f_raw(x(:).');
end
end

%% ========================================================================
% CHECKPOINT HELPERS
% ========================================================================
function LocalTaskTable = build_local_task_table(numAlgs, numR, linIdxMissing)
LocalTaskTable = zeros(numel(linIdxMissing), 3);
for k = 1:numel(linIdxMissing)
    [a, r] = ind2sub([numAlgs, numR], linIdxMissing(k));
    LocalTaskTable(k,:) = [a, r, linIdxMissing(k)];
end
end

function [RawFinal, TimeAlg, FEAlg, CurvesRaw, BestPos, DoneMask, OKMask, MsgMask, SeedMask] = ...
    merge_case_checkpoint(caseData, ci, RawFinal, TimeAlg, FEAlg, CurvesRaw, BestPos, DoneMask, OKMask, MsgMask, SeedMask)
[numAlg, Runs] = size(caseData.raw_ci);
for a = 1:numAlg
    for r = 1:Runs
        RawFinal(a,ci,r)   = caseData.raw_ci(a,r);
        TimeAlg(a,ci,r)    = caseData.time_ci(a,r);
        FEAlg(a,ci,r)      = caseData.fe_ci(a,r);
        CurvesRaw{a,ci,r}  = caseData.curves_ci{a,r};
        BestPos{a,ci,r}    = caseData.best_ci{a,r};
        OKMask(a,ci,r)     = caseData.ok_ci(a,r);
        MsgMask(a,ci,r)    = string(caseData.msg_ci(a,r));
        SeedMask(a,ci,r)   = caseData.seed_ci(a,r);
        DoneMask(a,ci,r)   = caseData.done_ci(a,r);
    end
end
end

function taskFile = build_task_checkpoint_file(taskDir, case_id, runIdx, algName)
safeAlg = regexprep(algName, '[^a-zA-Z0-9_]', '_');
taskFile = fullfile(taskDir, sprintf('C%02d_R%02d_%s.mat', case_id, runIdx, safeAlg));
end

function TaskResults = rebuild_task_results_from_arrays(RawFinal, TimeAlg, FEAlg, SeedMask, OKMask, MsgMask, Entries, Cases)
numAlgs = size(RawFinal,1); numCases = size(RawFinal,2); numRuns = size(RawFinal,3);
numTasks = numAlgs * numCases * numRuns;
TaskResults = table('Size', [numTasks 15], ...
    'VariableTypes', {'double','double','double','double','string','string','string','string','string','double','logical','double','double','double','string'}, ...
    'VariableNames', {'AlgIdx','CasePos','CaseID','RunIdx','Alias','Algorithm','Family','Role','PairID','Seed','OK','BestFit','TimeSec','FEs','Message'});
row = 0;
for a = 1:numAlgs
    for ci = 1:numCases
        for r = 1:numRuns
            row = row + 1;
            TaskResults.AlgIdx(row)    = a;
            TaskResults.CasePos(row)   = ci;
            TaskResults.CaseID(row)    = Cases(ci).id;
            TaskResults.RunIdx(row)    = r;
            TaskResults.Alias(row)     = string(Entries(a).alias);
            TaskResults.Algorithm(row) = string(Entries(a).name);
            TaskResults.Family(row)    = string(Entries(a).family);
            TaskResults.Role(row)      = string(Entries(a).role);
            TaskResults.PairID(row)    = string(Entries(a).pair_id);
            TaskResults.Seed(row)      = SeedMask(a,ci,r);
            TaskResults.OK(row)        = OKMask(a,ci,r);
            TaskResults.BestFit(row)   = RawFinal(a,ci,r);
            TaskResults.TimeSec(row)   = TimeAlg(a,ci,r);
            TaskResults.FEs(row)       = FEAlg(a,ci,r);
            TaskResults.Message(row)   = string(MsgMask(a,ci,r));
        end
    end
end
end

%% ========================================================================
% STATISTICS
% ========================================================================
function Summary = build_summary_table(Entries, RawFinal, TimeAlg, FEAlg, OKMask, DoneMask, RanksRunBlock)
numAlgs = numel(Entries);
Summary = table();
Summary.Alias             = strings(numAlgs,1);
Summary.Algorithm         = strings(numAlgs,1);
Summary.Family            = strings(numAlgs,1);
Summary.Role              = strings(numAlgs,1);
Summary.GlobalMeanFit     = nan(numAlgs,1);
Summary.GlobalStdFit      = nan(numAlgs,1);
Summary.GlobalMedianFit   = nan(numAlgs,1);
Summary.GlobalMeanTimeSec = nan(numAlgs,1);
Summary.GlobalMeanFEs     = nan(numAlgs,1);
Summary.ValidPct          = nan(numAlgs,1);
Summary.CompletedPct      = nan(numAlgs,1);
Summary.MeanRankRunCase   = nan(numAlgs,1);
Summary.FriedmanMeanRank  = nan(numAlgs,1);
for a = 1:numAlgs
    Summary.Alias(a)             = string(Entries(a).alias);
    Summary.Algorithm(a)         = string(Entries(a).name);
    Summary.Family(a)            = string(Entries(a).family);
    Summary.Role(a)              = string(Entries(a).role);
    Summary.GlobalMeanFit(a)     = mean(reshape(RawFinal(a,:,:),1,[]), 'omitnan');
    Summary.GlobalStdFit(a)      = std(reshape(RawFinal(a,:,:),1,[]), 0, 'omitnan');
    Summary.GlobalMedianFit(a)   = median(reshape(RawFinal(a,:,:),1,[]), 'omitnan');
    Summary.GlobalMeanTimeSec(a) = mean(reshape(TimeAlg(a,:,:),1,[]), 'omitnan');
    Summary.GlobalMeanFEs(a)     = mean(reshape(FEAlg(a,:,:),1,[]), 'omitnan');
    Summary.ValidPct(a)          = 100 * mean(reshape(OKMask(a,:,:),1,[]), 'omitnan');
    Summary.CompletedPct(a)      = 100 * mean(reshape(DoneMask(a,:,:),1,[]), 'omitnan');
    Summary.MeanRankRunCase(a)   = mean(reshape(RanksRunBlock(a,:,:),1,[]), 'omitnan');
end
end

function [Y_valid, p_friedman, MeanRanksFromFriedman] = build_friedman_from_fits(RawFinal)
[numAlgs, numCases, numRuns] = size(RawFinal);
Y = nan(numCases*numRuns, numAlgs);
row = 0;
for ci = 1:numCases
    for r = 1:numRuns
        row = row + 1;
        Y(row,:) = RawFinal(:,ci,r)';
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

function [MeanByCase, StdByCase, RankByMeanCase, AvgRank] = build_casewise_stats_from_mean(RawFinal)
[numAlgs, numCases, ~] = size(RawFinal);
MeanByCase = nan(numCases, numAlgs);
StdByCase = nan(numCases, numAlgs);
RankByMeanCase = nan(numCases, numAlgs);
for ci = 1:numCases
    for a = 1:numAlgs
        x = squeeze(RawFinal(a,ci,:));
        MeanByCase(ci,a) = mean(x, 'omitnan');
        StdByCase(ci,a)  = std(x, 0, 'omitnan');
    end
    RankByMeanCase(ci,:) = tiedrank_local(MeanByCase(ci,:)')';
end
AvgRank = mean(RankByMeanCase, 1, 'omitnan')';
end

function anchorIdx = pick_anchor_index_sarsa_by_meanrank_fir(Entries, AvgRankAll, Summary)
roles = string({Entries.role})';
cand = find(roles == "SARSA");
[~, ord] = sort(AvgRankAll(cand), 'ascend');
anchorIdx = cand(ord(1));

bestVal = AvgRankAll(anchorIdx);
tie = cand(abs(AvgRankAll(cand) - bestVal) <= 1e-12);
if numel(tie) > 1
    sumNames = string(Summary.Algorithm);
    vals = nan(numel(tie),1);
    for i = 1:numel(tie)
        idx = find(sumNames == string(Entries(tie(i)).name), 1, 'first');
        vals(i) = Summary.GlobalMeanFit(idx);
    end
    [~, j] = min(vals);
    anchorIdx = tie(j);
end
end

function Pairwise = build_pairwise_vs_anchor(Entries, RawFinal, anchorIdx)
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
allAnchor = reshape(RawFinal(anchorIdx,:,:), 1, []);
keep = 0; rawp = []; rowmap = [];
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    keep = keep + 1;
    x = reshape(RawFinal(a,:,:), 1, []);
    [p, ~, ~] = signrank_safe(x, allAnchor);
    Pairwise.Alias(keep) = string(Entries(a).alias);
    Pairwise.Algorithm(keep) = string(Entries(a).name);
    Pairwise.Family(keep) = string(Entries(a).family);
    Pairwise.Role(keep) = string(Entries(a).role);
    Pairwise.Anchor(keep) = string(Entries(anchorIdx).name);
    Pairwise.RawP(keep) = p;
    Pairwise.MedianDelta(keep) = median(x - allAnchor, 'omitnan');
    Pairwise.BetterThanAnchorPct(keep) = 100 * mean(x < allAnchor, 'omitnan');
    rawp(end+1,1) = p; %#ok<AGROW>
    rowmap(end+1,1) = keep; %#ok<AGROW>
end
Pairwise = Pairwise(1:keep,:);
validP = ~isnan(rawp);
adjp_all = nan(size(rawp)); rej_all = false(size(rawp));
if any(validP)
    [adjp_tmp, rej_tmp] = holm_correction(rawp(validP), 0.05);
    adjp_all(validP) = adjp_tmp; rej_all(validP) = rej_tmp;
end
for i = 1:numel(rowmap)
    Pairwise.HolmAdjP(rowmap(i)) = adjp_all(i);
    Pairwise.RejectHolm(rowmap(i)) = rej_all(i);
end
Pairwise = sortrows(Pairwise, 'RawP', 'ascend');
end

function Tbl = build_wtl_by_case(RawFinal, Entries, anchorIdx, Cases, tieTol)
numAlgs = size(RawFinal,1); numCases = size(RawFinal,2); maxRows = (numAlgs-1)*numCases;
Tbl = table('Size', [maxRows 7], ...
    'VariableTypes', {'string','string','string','string','double','double','string'}, ...
    'VariableNames', {'Alias','Algorithm','Family','Role','CaseID','AnchorIndex','WTL'});
row = 0;
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    for ci = 1:numCases
        row = row + 1;
        x = squeeze(RawFinal(a,ci,:));
        y = squeeze(RawFinal(anchorIdx,ci,:));
        [w,t,l] = count_wtl(x,y,tieTol);
        Tbl.Alias(row) = string(Entries(a).alias);
        Tbl.Algorithm(row) = string(Entries(a).name);
        Tbl.Family(row) = string(Entries(a).family);
        Tbl.Role(row) = string(Entries(a).role);
        Tbl.CaseID(row) = Cases(ci).id;
        Tbl.AnchorIndex(row) = anchorIdx;
        Tbl.WTL(row) = sprintf('%d/%d/%d', w,t,l);
    end
end
Tbl = Tbl(1:row,:);
end

function Tbl = build_wtl_global(RawFinal, Entries, anchorIdx, tieTol)
numAlgs = size(RawFinal,1); rows = numAlgs - 1;
Tbl = table('Size', [rows 6], ...
    'VariableTypes', {'string','string','string','string','string','double'}, ...
    'VariableNames', {'Alias','Algorithm','Family','Role','WTL_Global','AnchorIndex'});
row = 0; anchor = reshape(RawFinal(anchorIdx,:,:), 1, []);
for a = 1:numAlgs
    if a == anchorIdx, continue; end
    row = row + 1;
    x = reshape(RawFinal(a,:,:),1,[]);
    [w,t,l] = count_wtl(x, anchor, tieTol);
    Tbl.Alias(row) = string(Entries(a).alias);
    Tbl.Algorithm(row) = string(Entries(a).name);
    Tbl.Family(row) = string(Entries(a).family);
    Tbl.Role(row) = string(Entries(a).role);
    Tbl.WTL_Global(row) = sprintf('%d/%d/%d', w,t,l);
    Tbl.AnchorIndex(row) = anchorIdx;
end
Tbl = Tbl(1:row,:);
end

function [w,t,l] = count_wtl(x,y,tieTol)
mask = ~(isnan(x) | isnan(y));
x = x(mask); y = y(mask);
if isempty(x), w=0; t=0; l=0; return; end
d = x - y;
w = sum(d < -tieTol);
t = sum(abs(d) <= tieTol);
l = sum(d > tieTol);
end

%% ========================================================================
% EXPORT MAIN-TEXT TABLES (UNIFIED)
% ========================================================================
function export_fir_ablation_tables(Entries, Pairwise, WTL_global, Cases, ...
    anchorIdx, anchorAlias, MeanByCaseAll, StdByCaseAll, RankByMeanCaseAll, AvgRankAll, tables_dir, csv_dir)
idxKeep = 1:numel(Entries);
ItemLabels = "C" + string([Cases.id]);
AliasTbl = build_alias_table_ablation(Entries);
MainTbl = build_sensitivity_panel_table_unified(Entries, idxKeep, ItemLabels, 'Case', MeanByCaseAll, StdByCaseAll, RankByMeanCaseAll, AvgRankAll);
StatsTbl = build_anchor_stats_table_ablation(Entries, Pairwise, WTL_global, AvgRankAll, anchorIdx);

writetable(AliasTbl, fullfile(csv_dir, 'fir_ablation_alias_map.csv'));
writetable(MainTbl,  fullfile(csv_dir, 'fir_ablation_main.csv'));
writetable(StatsTbl, fullfile(csv_dir, 'fir_ablation_anchor_stats.csv'));

write_latex_alias_map_ablation(AliasTbl, fullfile(tables_dir, 'fir_ablation_alias_map.tex'));
write_latex_ablation_panel_unified(Entries, idxKeep, ItemLabels, 'Case', 'FIR', 'fitness', MeanByCaseAll, StdByCaseAll, RankByMeanCaseAll, AvgRankAll, ...
    fullfile(tables_dir, 'fir_ablation_main.tex'));
write_latex_anchor_stats_unified(StatsTbl, string(anchorAlias), fullfile(tables_dir, 'fir_ablation_anchor_stats.tex'));
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
        case 'MPHB'
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
fprintf(fid, '\\caption{Alias map for the FIR ablation study.}\n');
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
function export_convergence_panels_fir(CurvesRaw, FEAlg, Entries, Cases, fig_dir)
mkdir_if_needed(fig_dir);

labels = cell(numel(Entries),1);
for a = 1:numel(Entries)
    labels{a} = clean_plot_label(Entries(a).alias);
end

[~, ord] = sort(cellfun(@plot_sort_key, labels));
labelsOrd = labels(ord);

fig = figure('Visible','off','Color','w','Units','pixels','Position',[50 50 1520 800]);
t = tiledlayout(ceil(numel(Cases)/4), 4, 'TileSpacing','tight', 'Padding','compact');
legendHandles = gobjects(numel(ord),1);

for ci = 1:numel(Cases)
    ax = nexttile(t);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    for kk = 1:numel(ord)
        a = ord(kk);
        [xfe, ymean] = build_mean_curve_from_cells(CurvesRaw(a,ci,:), FEAlg(a,ci,:));
        if isempty(xfe), continue; end
        yplot = max(ymean(:).', eps);
        [mk, ls] = get_plot_style(kk);
        idxMarks = unique(round(linspace(1, numel(xfe), min(8, numel(xfe)))));
        h = plot(ax, xfe, yplot, 'LineWidth', 1.0, ...
            'LineStyle', ls, 'Marker', mk, 'MarkerIndices', idxMarks, ...
            'MarkerSize', 5);
        if ci == 1
            legendHandles(kk) = h;
        end
    end

    set(ax, 'XScale','linear', 'YScale','log');
    xlabel(ax, 'FEs', 'FontSize', 11);
    ylabel(ax, 'Mean fitness', 'FontSize', 11);
    title(ax, sprintf('C%d - %s', Cases(ci).id, Cases(ci).name), 'Interpreter','none', 'FontSize', 12, 'FontWeight', 'bold');
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

saveas(fig, fullfile(fig_dir, 'fir_convergence_panels_logFE.png'));
savefig(fig, fullfile(fig_dir, 'fir_convergence_panels_logFE.fig'));
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
fprintf(fid, '\\caption{Alias map used throughout the FIR sensitivity and ablation analysis.}\n');
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


function write_entries_manifest(Entries, outcsv)
T = table();
T.Alias     = string({Entries.alias})';
T.PairID    = string({Entries.pair_id})';
T.CFG_ID    = string({Entries.cfg_id})';
T.Algorithm = string({Entries.name})';
T.Family    = string({Entries.family})';
T.Role      = string({Entries.role})';
T.Label     = string({Entries.label})';

n = numel(Entries);
Kvec   = nan(n,1);
Nivec  = nan(n,1);
Gvec   = nan(n,1);
Wvec   = nan(n,1);
Rvec   = nan(n,1);
Tvec   = nan(n,1);
FEvec  = nan(n,1);
Bvec   = false(n,1);
Avec   = false(n,1);

for i = 1:n
    Kvec(i)  = Entries(i).K;
    Nivec(i) = Entries(i).Ni;
    Gvec(i)  = Entries(i).gamma;
    Wvec(i)  = Entries(i).w_best;
    Rvec(i)  = Entries(i).rho_sub;
    Tvec(i)  = Entries(i).T;
    FEvec(i) = Entries(i).EstimatedFEs;
    Bvec(i)  = logical(Entries(i).fad_before_p2);
    Avec(i)  = logical(Entries(i).fad_after_p2);
end

T.K             = Kvec;
T.Ni            = Nivec;
T.gamma         = Gvec;
T.w_best        = Wvec;
T.rho_sub       = Rvec;
T.fad_before_p2 = Bvec;
T.fad_after_p2  = Avec;
T.T             = Tvec;
T.EstimatedFEs  = FEvec;

writetable(T, outcsv);
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
