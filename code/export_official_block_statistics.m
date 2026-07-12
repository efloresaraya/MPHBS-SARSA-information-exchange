function export_official_block_statistics()
%EXPORT_OFFICIAL_BLOCK_STATISTICS Reproduce the manuscript statistics.
%
% Each CEC2022 function or FIR case is treated as one independent
% benchmark block. The 30 common-seed runs remain paired within each
% block. The export includes:
%   - paired Wilcoxon signed-rank tests within each function/case;
%   - domain-wide Holm correction over all anchor-comparator instances;
%   - Friedman omnibus tests over function/case mean ranks;
%   - Kendall's W as the omnibus effect size; and
%   - anchor-based Friedman post hoc comparisons with Holm correction.

codeDir = fileparts(mfilename('fullpath'));
packageDir = fileparts(codeDir);
resultsDir = fullfile(packageDir, 'results');

campaigns = {
    'CEC2022_sensitivity_common_seeds', ...
        'CEC2022_SENSITIVITY', 'B-C1', 'FunctionID', 'Error';
    'CEC2022_ablation_BC1_common_seeds', ...
        'CEC2022_ABLATION', 'BS-C1', 'FunctionID', 'Error';
    'CEC2022_external_BC1_common_seeds', ...
        'CEC2022_EXTERNAL', 'B-C1', 'FunctionID', 'Error';
    'FIR_sensitivity_common_seeds', ...
        'FIR_SENSITIVITY', 'A-C1', 'CaseID', 'BestFit';
    'FIR_ablation_AC1_common_seeds', ...
        'FIR_ABLATION', 'AS-C1', 'CaseID', 'BestFit';
    'FIR_external_AC1_common_seeds', ...
        'FIR_EXTERNAL', 'A-C1', 'CaseID', 'BestFit'
};

summaryRows = cell(size(campaigns, 1), 1);
for campaignIdx = 1:size(campaigns, 1)
    resultDir = fullfile(resultsDir, campaigns{campaignIdx, 1});
    csvFile = fullfile(resultDir, 'csv', 'task_results.csv');
    outputDir = fullfile(resultDir, 'statistics_official_block_level');
    summaryRows{campaignIdx} = analyze_campaign( ...
        csvFile, campaigns{campaignIdx, 2}, campaigns{campaignIdx, 3}, ...
        campaigns{campaignIdx, 4}, campaigns{campaignIdx, 5}, outputDir);
    remove_runner_native_inference(resultDir);
end

masterSummary = vertcat(summaryRows{:});
writetable(masterSummary, fullfile(packageDir, 'documentation', ...
    'OFFICIAL_OMNIBUS_STATISTICS.csv'));

fprintf('\nOfficial block-level statistics exported for %d campaigns.\n', ...
    size(campaigns, 1));
end

function remove_runner_native_inference(resultDir)
% Keep the public output focused on the official block-level analysis.
csvFiles = [dir(fullfile(resultDir, 'csv', 'pairwise_vs_anchor.csv')); ...
    dir(fullfile(resultDir, 'csv', '*_anchor_stats.csv'))];
for fileIdx = 1:numel(csvFiles)
    delete(fullfile(csvFiles(fileIdx).folder, csvFiles(fileIdx).name));
end

tableFiles = dir(fullfile(resultDir, 'tables', '*_anchor_stats.tex'));
for fileIdx = 1:numel(tableFiles)
    delete(fullfile(tableFiles(fileIdx).folder, tableFiles(fileIdx).name));
end

oldDirectory = fullfile(resultDir, 'statistics_paired_complete');
if isfolder(oldDirectory)
    rmdir(oldDirectory, 's');
end
end

function omnibus = analyze_campaign(csvFile, domainName, anchorAlias, ...
        instanceVariable, valueVariable, outputDir)
if ~isfile(csvFile)
    error('Missing campaign data: %s', csvFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

data = readtable(csvFile, 'TextType', 'string');
data = data(data.OK == 1, :);
aliases = unique(data.Alias, 'stable');
competitors = aliases(aliases ~= anchorAlias);
instances = unique(data.(instanceVariable), 'sorted');

assert(any(aliases == anchorAlias), 'Anchor %s was not found.', anchorAlias);
assert(all(isfinite(data.(valueVariable))), ...
    '%s contains nonfinite objective values.', domainName);

instanceTests = paired_instance_tests(data, domainName, anchorAlias, ...
    competitors, instances, instanceVariable, valueVariable);
instanceTests.HolmDomainP = holm_adjust(instanceTests.RawP);
instanceTests.RejectHolmDomain = instanceTests.HolmDomainP < 0.05;
writetable(instanceTests, fullfile(outputDir, ...
    'official_paired_wilcoxon_by_instance.csv'));

[meanMatrix, rankMatrix] = build_instance_matrices(data, aliases, ...
    instances, instanceVariable, valueVariable);
[friedmanP, friedmanTable, friedmanStats] = friedman( ...
    meanMatrix, 1, 'off');
avgRanks = mean(rankMatrix, 1);
kendallW = kendall_concordance(rankMatrix);

anchorIndex = find(aliases == anchorAlias, 1);
globalPairwise = block_level_posthoc(data, aliases, competitors, ...
    anchorAlias, anchorIndex, instances, meanMatrix, avgRanks, ...
    instanceTests);
writetable(globalPairwise, fullfile(outputDir, ...
    'official_global_by_instance_blocks.csv'));

rankTable = array2table(rankMatrix, ...
    'VariableNames', matlab.lang.makeValidName(cellstr(aliases)));
rankTable = addvars(rankTable, instances, 'Before', 1, ...
    'NewVariableNames', 'Instance');
writetable(rankTable, fullfile(outputDir, ...
    'official_instance_ranks.csv'));

Domain = string(domainName);
Anchor = string(anchorAlias);
NInstances = numel(instances);
NMethods = numel(aliases);
FriedmanChiSquare = friedmanTable{2, 5};
FriedmanP = friedmanP;
KendallW = kendallW;
AnchorAvgRank = avgRanks(anchorIndex);
omnibus = table(Domain, Anchor, NInstances, NMethods, FriedmanChiSquare, ...
    FriedmanP, KendallW, AnchorAvgRank);
writetable(omnibus, fullfile(outputDir, ...
    'official_friedman_kendall_omnibus.csv'));
write_latex_global_table(globalPairwise, omnibus, outputDir);

save(fullfile(outputDir, 'official_block_statistics.mat'), ...
    'domainName', 'anchorAlias', 'aliases', 'competitors', 'instances', ...
    'meanMatrix', 'rankMatrix', 'avgRanks', 'friedmanP', ...
    'friedmanTable', 'friedmanStats', 'kendallW', 'instanceTests', ...
    'globalPairwise', 'omnibus');

fprintf('%-22s Friedman p=%-12.6g Kendall W=%.6f anchor rank=%.4f\n', ...
    domainName, friedmanP, kendallW, avgRanks(anchorIndex));
end

function write_latex_global_table(globalPairwise, omnibus, outputDir)
fileName = fullfile(outputDir, ...
    'official_global_by_instance_blocks.tex');
fileId = fopen(fileName, 'w');
assert(fileId >= 0, 'Could not write %s.', fileName);
cleanup = onCleanup(@() fclose(fileId));

domain = latex_escape(globalPairwise.Domain(1));
anchor = latex_escape(globalPairwise.Anchor(1));
fprintf(fileId, '\\begin{table*}[!t]\n');
fprintf(fileId, '\\centering\n');
fprintf(fileId, '\\footnotesize\n');
fprintf(fileId, ['\\caption{%s comparison using benchmark instances as ' ...
    'independent blocks (anchor = %s). Friedman $p=%.3e$ and ' ...
    'Kendall''s $W=%.3f$. Post hoc $p$-values are Holm-adjusted ' ...
    'across competitors.}\n'], domain, anchor, omnibus.FriedmanP, ...
    omnibus.KendallW);
fprintf(fileId, '\\begin{tabular}{llcccccc}\n');
fprintf(fileId, '\\hline\n');
fprintf(fileId, ['Method & Role & W/T/L & Anchor rank & Comp. rank & ' ...
    'Holm $p$ & Reject & Sig. inst. A/C ']);
fprintf(fileId, '%s\n', '\\');
fprintf(fileId, '\\hline\n');
for rowIdx = 1:height(globalPairwise)
    row = globalPairwise(rowIdx, :);
    method = latex_escape(row.Competitor);
    role = latex_escape(row.Role);
    wtl = sprintf('%d/%d/%d', row.FunctionWinsAnchor, ...
        row.FunctionTies, row.FunctionLossesAnchor);
    significant = sprintf('%d/%d', row.SignificantInstancesAnchor, ...
        row.SignificantInstancesCompetitor);
    fprintf(fileId, ...
        '%s & %s & %s & %.3f & %.3f & %.3e & %d & %s ', ...
        method, role, wtl, row.AnchorAvgRank, row.CompetitorAvgRank, ...
        row.FriedmanPostHocHolmP, row.RejectFriedmanHolm, significant);
    fprintf(fileId, '%s\n', '\\');
end
fprintf(fileId, '\\hline\n');
fprintf(fileId, '\\end{tabular}\n');
fprintf(fileId, '\\end{table*}\n');
clear cleanup
end

function escaped = latex_escape(value)
escaped = char(string(value));
escaped = strrep(escaped, '\', '\textbackslash{}');
escaped = strrep(escaped, '_', '\_');
escaped = strrep(escaped, '&', '\&');
escaped = strrep(escaped, '%', '\%');
escaped = strrep(escaped, '#', '\#');
end

function tests = paired_instance_tests(data, domainName, anchorAlias, ...
        competitors, instances, instanceVariable, valueVariable)
nTests = numel(competitors) * numel(instances);
Domain = repmat(string(domainName), nTests, 1);
Instance = zeros(nTests, 1);
Anchor = repmat(string(anchorAlias), nTests, 1);
Competitor = strings(nTests, 1);
NPairs = zeros(nTests, 1);
WinsAnchor = zeros(nTests, 1);
Ties = zeros(nTests, 1);
LossesAnchor = zeros(nTests, 1);
MeanAnchor = zeros(nTests, 1);
MeanCompetitor = zeros(nTests, 1);
MedianAnchor = zeros(nTests, 1);
MedianCompetitor = zeros(nTests, 1);
RawP = ones(nTests, 1);

selectedVariables = ["RunIdx", "Seed", string(valueVariable)];
row = 0;
for competitorIdx = 1:numel(competitors)
    competitor = competitors(competitorIdx);
    for instanceIdx = 1:numel(instances)
        row = row + 1;
        instance = instances(instanceIdx);
        anchorData = data(data.Alias == anchorAlias & ...
            data.(instanceVariable) == instance, selectedVariables);
        competitorData = data(data.Alias == competitor & ...
            data.(instanceVariable) == instance, selectedVariables);
        anchorData.Properties.VariableNames{3} = 'AnchorValue';
        competitorData.Properties.VariableNames{3} = 'CompetitorValue';
        paired = innerjoin(anchorData, competitorData, ...
            'Keys', {'RunIdx', 'Seed'});

        assert(height(paired) == 30, ...
            '%s %s instance %d has %d paired observations.', ...
            domainName, competitor, instance, height(paired));
        assert(numel(unique(paired.Seed)) == 30, ...
            '%s %s instance %d does not contain 30 unique seeds.', ...
            domainName, competitor, instance);

        anchorValues = paired.AnchorValue;
        competitorValues = paired.CompetitorValue;
        tolerance = comparison_tolerance(anchorValues, competitorValues);
        difference = competitorValues - anchorValues;
        testDifference = difference;
        testDifference(abs(testDifference) <= tolerance) = 0;

        Instance(row) = instance;
        Competitor(row) = competitor;
        NPairs(row) = height(paired);
        WinsAnchor(row) = sum(difference > tolerance);
        Ties(row) = sum(abs(difference) <= tolerance);
        LossesAnchor(row) = sum(difference < -tolerance);
        MeanAnchor(row) = mean(anchorValues);
        MeanCompetitor(row) = mean(competitorValues);
        MedianAnchor(row) = median(anchorValues);
        MedianCompetitor(row) = median(competitorValues);
        if ~all(testDifference == 0)
            RawP(row) = signrank(testDifference, 0, ...
                'method', 'approximate');
        end
    end
end

tests = table(Domain, Instance, Anchor, Competitor, NPairs, WinsAnchor, ...
    Ties, LossesAnchor, MeanAnchor, MeanCompetitor, MedianAnchor, ...
    MedianCompetitor, RawP);
end

function [meanMatrix, rankMatrix] = build_instance_matrices(data, aliases, ...
        instances, instanceVariable, valueVariable)
meanMatrix = zeros(numel(instances), numel(aliases));
rankMatrix = zeros(size(meanMatrix));
for instanceIdx = 1:numel(instances)
    for aliasIdx = 1:numel(aliases)
        values = data.(valueVariable)( ...
            data.(instanceVariable) == instances(instanceIdx) & ...
            data.Alias == aliases(aliasIdx));
        assert(numel(values) == 30, ...
            'Instance %d, method %s has %d runs instead of 30.', ...
            instances(instanceIdx), aliases(aliasIdx), numel(values));
        meanMatrix(instanceIdx, aliasIdx) = mean(values);
    end
    rankMatrix(instanceIdx, :) = tiedrank(meanMatrix(instanceIdx, :));
end
end

function globalPairwise = block_level_posthoc(data, aliases, competitors, ...
        anchorAlias, anchorIndex, instances, meanMatrix, avgRanks, ...
        instanceTests)
nCompetitors = numel(competitors);
nInstances = numel(instances);
nMethods = numel(aliases);
standardError = sqrt(nMethods * (nMethods + 1) / (6 * nInstances));

Domain = repmat(instanceTests.Domain(1), nCompetitors, 1);
Anchor = repmat(string(anchorAlias), nCompetitors, 1);
Competitor = competitors;
Role = strings(nCompetitors, 1);
FunctionWinsAnchor = zeros(nCompetitors, 1);
FunctionTies = zeros(nCompetitors, 1);
FunctionLossesAnchor = zeros(nCompetitors, 1);
AnchorAvgRank = repmat(avgRanks(anchorIndex), nCompetitors, 1);
CompetitorAvgRank = zeros(nCompetitors, 1);
FriedmanPostHocRawP = zeros(nCompetitors, 1);
SignificantInstancesAnchor = zeros(nCompetitors, 1);
SignificantInstancesCompetitor = zeros(nCompetitors, 1);

for competitorIdx = 1:nCompetitors
    competitor = competitors(competitorIdx);
    methodRows = data(data.Alias == competitor, :);
    Role(competitorIdx) = methodRows.Role(1);
    methodIndex = find(aliases == competitor, 1);
    delta = meanMatrix(:, methodIndex) - meanMatrix(:, anchorIndex);
    tolerance = comparison_tolerance( ...
        meanMatrix(:, anchorIndex), meanMatrix(:, methodIndex));

    FunctionWinsAnchor(competitorIdx) = sum(delta > tolerance);
    FunctionTies(competitorIdx) = sum(abs(delta) <= tolerance);
    FunctionLossesAnchor(competitorIdx) = sum(delta < -tolerance);
    CompetitorAvgRank(competitorIdx) = avgRanks(methodIndex);

    zValue = (avgRanks(anchorIndex) - avgRanks(methodIndex)) / ...
        standardError;
    FriedmanPostHocRawP(competitorIdx) = ...
        2 * normcdf(-abs(zValue));

    rows = instanceTests.Competitor == competitor;
    significantRows = rows & instanceTests.RejectHolmDomain;
    SignificantInstancesAnchor(competitorIdx) = sum(significantRows & ...
        instanceTests.MeanAnchor < instanceTests.MeanCompetitor);
    SignificantInstancesCompetitor(competitorIdx) = sum(significantRows & ...
        instanceTests.MeanCompetitor < instanceTests.MeanAnchor);
end

FriedmanPostHocHolmP = holm_adjust(FriedmanPostHocRawP);
RejectFriedmanHolm = FriedmanPostHocHolmP < 0.05;
globalPairwise = table(Domain, Anchor, Competitor, Role, ...
    FunctionWinsAnchor, FunctionTies, FunctionLossesAnchor, ...
    AnchorAvgRank, CompetitorAvgRank, FriedmanPostHocRawP, ...
    FriedmanPostHocHolmP, RejectFriedmanHolm, ...
    SignificantInstancesAnchor, SignificantInstancesCompetitor);
end

function tolerance = comparison_tolerance(firstValues, secondValues)
tolerance = max(1e-14, 1e-12 * ...
    max([abs(firstValues(:)); abs(secondValues(:)); 1]));
end

function w = kendall_concordance(rankMatrix)
[nBlocks, nMethods] = size(rankMatrix);
rankSums = sum(rankMatrix, 1);
center = nBlocks * (nMethods + 1) / 2;
sumSquares = sum((rankSums - center).^2);

tieCorrection = 0;
for blockIdx = 1:nBlocks
    [~, ~, groups] = unique(rankMatrix(blockIdx, :));
    groupSizes = accumarray(groups(:), 1);
    tieCorrection = tieCorrection + ...
        sum(groupSizes.^3 - groupSizes);
end

denominator = nBlocks^2 * (nMethods^3 - nMethods) - ...
    nBlocks * tieCorrection;
if denominator == 0
    w = 0;
else
    w = 12 * sumSquares / denominator;
end
end

function adjusted = holm_adjust(pValues)
pValues = pValues(:);
nValues = numel(pValues);
[sortedP, order] = sort(pValues, 'ascend');
sortedAdjusted = zeros(nValues, 1);
runningMaximum = 0;
for index = 1:nValues
    runningMaximum = max(runningMaximum, ...
        (nValues - index + 1) * sortedP(index));
    sortedAdjusted(index) = min(runningMaximum, 1);
end
adjusted = zeros(nValues, 1);
adjusted(order) = sortedAdjusted;
end
