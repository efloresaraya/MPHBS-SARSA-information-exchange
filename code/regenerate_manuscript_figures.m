% Regenerate the four manuscript convergence figures from final paired data.
% Only presentation properties and paper-facing aliases are changed.

packageDir = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(packageDir, 'results');
manuscriptDir = fullfile(packageDir, 'manuscript_figures');

jobs = struct([]);
jobs(1).source = fullfile(resultsDir, ...
    'CEC2022_sensitivity_common_seeds', ...
    'figures', 'cec_convergence_panels_logFE.fig');
jobs(1).output = fullfile(manuscriptDir, 'cec', ...
    'cec_internal_sensitivity_convergence_logFE');
jobs(1).domain = 'cec';
jobs(1).comparison = 'internal';

jobs(2).source = fullfile(resultsDir, ...
    'FIR_sensitivity_common_seeds', ...
    'figures', 'fir_convergence_panels_logFE.fig');
jobs(2).output = fullfile(manuscriptDir, 'fir', ...
    'fir_internal_sensitivity_convergence_logFE');
jobs(2).domain = 'fir';
jobs(2).comparison = 'internal';

jobs(3).source = fullfile(resultsDir, ...
    'CEC2022_external_BC1_common_seeds', ...
    'figures', 'cec_convergence_panels_logFE.fig');
jobs(3).output = fullfile(manuscriptDir, 'cec', ...
    'cec_external_convergence_logFE');
jobs(3).domain = 'cec';
jobs(3).comparison = 'external';

jobs(4).source = fullfile(resultsDir, ...
    'FIR_external_AC1_common_seeds', ...
    'figures', 'fir_convergence_panels_logFE.fig');
jobs(4).output = fullfile(manuscriptDir, 'fir', ...
    'fir_external_convergence_logFE');
jobs(4).domain = 'fir';
jobs(4).comparison = 'external';

for jobIndex = 1:numel(jobs)
    job = jobs(jobIndex);
    assert(isfile(job.source), 'Missing source figure: %s', job.source);

    fig = openfig(job.source, 'invisible');
    cleanup = onCleanup(@() close_if_valid(fig));
    set(fig, 'Color', 'white');
    visualScale = 1;
    if strcmp(job.domain, 'fir')
        visualScale = 1.25;
        set(fig, 'Units', 'pixels');
        figurePosition = fig.Position;
        figurePosition(3:4) = [1280, 740];
        fig.Position = figurePosition;
    end

    axesHandles = findall(fig, 'Type', 'axes');
    for axisIndex = 1:numel(axesHandles)
        ax = axesHandles(axisIndex);
        set(ax, 'FontSize', 12 * visualScale, ...
            'LineWidth', 0.9 * visualScale);
        ax.XLabel.FontSize = 13 * visualScale;
        ax.YLabel.FontSize = 13 * visualScale;
        ax.Title.FontSize = 14 * visualScale;
        ax.Title.FontWeight = 'bold';

        if strcmp(job.domain, 'fir')
            ax.PlotBoxAspectRatio = [1, 1, 1];
            ax.PlotBoxAspectRatioMode = 'manual';
            titleText = string(ax.Title.String);
            titleText = regexprep(titleText, '^C(\d+)\s*-\s*', 'FIR$1 - ');
            ax.Title.String = char(titleText);
        end

        lineHandles = findall(ax, 'Type', 'line');
        for lineIndex = 1:numel(lineHandles)
            lineHandles(lineIndex).LineWidth = max( ...
                1.5 * visualScale, lineHandles(lineIndex).LineWidth);
            lineHandles(lineIndex).MarkerSize = max( ...
                6 * visualScale, lineHandles(lineIndex).MarkerSize);
        end
    end

    legendHandles = findall(fig, 'Type', 'legend');
    assert(numel(legendHandles) == 1, ...
        'Expected one legend in %s, found %d.', job.source, numel(legendHandles));
    lgd = legendHandles(1);
    labels = string(lgd.String);
    labels = paper_aliases(labels);
    lgd.String = cellstr(labels);
    lgd.Interpreter = 'none';
    lgd.Orientation = 'horizontal';

    if strcmp(job.comparison, 'external')
        lgd.NumColumns = 5;
        lgd.FontSize = 12 * visualScale;
    else
        lgd.NumColumns = 8;
        lgd.FontSize = 10.5 * visualScale;
    end

    try
        lgd.Layout.Tile = 'south';
    catch
        lgd.Location = 'southoutside';
    end

    drawnow;
    savefig(fig, job.output + ".fig");
    exportgraphics(fig, job.output + ".png", ...
        'Resolution', 400, 'BackgroundColor', 'white');
    fprintf('Generated %s.[fig|png]\n', job.output);
    clear cleanup
end

function labels = paper_aliases(labels)
    labels = replace(labels, "MPHBS-BC", "B-C");
    labels = replace(labels, "MPHBS-AC", "A-C");
    labels = replace(labels, "MPHBR-RBC", "RB-C");
    labels = replace(labels, "MPHBR-RAC", "RA-C");
    labels = replace(labels, "FDBTLABC", "FDB-TLABC");
end

function close_if_valid(fig)
    if isgraphics(fig)
        close(fig);
    end
end
