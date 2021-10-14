%% Load data

load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed/reference_bulk_profiles.mat', 'ref', 'r_bulk')
GM12878 = load('data/processed/GM12878.mat', 'replication_state_filtered', ...
    'percent_replicated_filtered', 'is_included_chr', 'aggregate_S_G1');

GM18507 = load('data/processed/GM18507.mat', 'replication_state_filtered', ...
    'percent_replicated_filtered', 'is_included_chr', 'aggregate_S_G1');

%% Figure skeleton
figure2 = figure;
set(figure2, 'Position', [25 9 6.5 8.2])

panelA = struct('top', axes('Units', 'inches', 'Position', [0.45 7.1278 5.8 0.5]), ...
'bottom', axes('Units', 'inches', 'Position', [0.45 5.5778 5.8 1.5]));

inset = axes('Units', 'inches', 'Position', [0.45 4.8 5.8 0.7778]);

panelA1 = struct('top', axes('Units', 'inches', 'Position', [0.75 4.3 0.8 0.5]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.75 3 0.8 1.25]));

panelA2 = struct('top', axes('Units', 'inches', 'Position', [2.75 4.3 3.5 0.5]), ...
    'bottom', axes('Units', 'inches', 'Position', [2.75 3 3.5 1.25]));

panelB = struct('top', axes('Units', 'inches', 'Position', [0.45 1.69 5.8 0.5]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.45 0.44 5.8 1.2]));

%% Heatmaps

Chr = 2;

params = struct('panel', {panelA panelA1 panelA2 panelB}, 'Chr', 2, 'X', [100 235], ...
    'sample', GM12878, 'title', '', 'inset', [NaN NaN], 'name', 'GM12878');

params(1).title = 'GM12878, 10x Genomics';
params(2).X = [125 130];
params(2).inset = [107.1 125.5];
params(3).X = [200 230];
params(3).inset = [153.5 235];
params(4).sample = GM18507;
params(4).title = 'GM18507, DLP+';
params(4).name = 'GM18507';

for p = 1:4
    
    top = params(p).panel.top;
    bottom = params(p).panel.bottom;
    
    plot(ref.(params(p).name){Chr}(:, 1) ./ 1e6, ref.(params(p).name){Chr}(:, 2), ...
        '.', 'Color', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', top)
    plot(params(p).sample.aggregate_S_G1{Chr}(:, 1) ./ 1e6, ...
        params(p).sample.aggregate_S_G1{Chr}(:, 2), '.', 'Color', s_light{1}, ...
        'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', top)
    set(top, 'XLim', params(p).X, 'XTick', [], 'YTick', [-1 1])
    ylabel(top, 'RT')
    title(top, params(p).title)
    
    if ismember(p, [1 4])
        yyaxis(top, 'right')
        set(top, 'YColor', 'k', 'YTick', [])
        ylabel(top, ['r = ' num2str(r_bulk.(params(p).name), '%0.2f')])
    end
    
    num_cells = sum(params(p).sample.is_included_chr(Chr, :));
    index = params(p).sample.percent_replicated_filtered(params(p).sample.is_included_chr(Chr, :));
    [Yticks, YLabels] = get_heatmap_yticks(index);
    
    r = params(p).sample.replication_state_filtered{Chr}(:, params(p).sample.is_included_chr(Chr, :));
    imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
        'Parent', bottom);
    set(bottom, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [1 num_cells], 'CLim', [2 4], ...
        'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    xlabel(bottom, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'], 'BackgroundColor', 'white')
    ylabel(bottom, '% Replicated')
    colormap(bottom, [convert_hex(g1_light); convert_hex(s_light{1})])
    yyaxis(bottom, 'right')
    set(bottom, 'YColor', 'k', 'YTick', [])
    ylabel(bottom, [num2str(num_cells) '  cells'])
    
end

xlabel(panelA1.bottom, ['Chr. ' num2str(Chr) ', Mb'], 'BackgroundColor', 'white')

%% Legend for panel A

plot(0, 0, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Bulk-seq', 'Parent', panelA.top)
plot(0, 0, 'Color', s_light{1}, 'LineWidth', 2, 'DisplayName', 'Single-Cell S/G1 Aggregate', ...
    'Parent', panelA.top)
legendA = legend(panelA.top);
set(legendA, 'Orientation', 'horizontal', 'FontSize', 9, 'Position', [0.0673 0.965 0.4241 0.02])
legendA.ItemTokenSize(1) = 15;

colorbarA = colorbar(panelA.bottom);
set(colorbarA, 'Orientation', 'horizontal', 'Position', [0.71 0.965 0.15 0.02],  ...
    'Ticks', [2.5 3.5], 'TickLabels', [2 4], 'YAxisLocation', 'bottom')

%% Connect insets

set(inset, 'XLim', params(1).X, 'YLim', [0 1], 'Visible', 'off')
uistack(inset, 'bottom')

for p = 2:3
    plot(params(p).X(1) * ones(1, 2), [-2 2], 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off', ...
        'Parent', panelA.top)
    plot(params(p).X(2) * ones(1, 2), [-2 2], 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off', ...
        'Parent', panelA.top)

    plot(params(p).X(1) * ones(1, 2), [1 1998], 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off', ...
        'Parent', panelA.bottom)
    plot(params(p).X(2) * ones(1, 2), [1 1998], 'k-', 'LineWidth', 0.7, 'HandleVisibility', 'off', ...
        'Parent', panelA.bottom)
    
    plot([params(p).X(1) params(p).inset(1)], [1 0], 'k-', 'LineWidth', 0.7, 'Parent', inset)
    plot([params(p).X(2) params(p).inset(2)], [1 0], 'k-', 'LineWidth', 0.7, 'Parent', inset)
end

%% Annotate panels
 
params = struct('panel', {panelA.top panelB.top}, 'text', {'a', 'b'});

for p = 1:2
    text(-0.0503, 1.4167, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/Figure2.pdf')
close
