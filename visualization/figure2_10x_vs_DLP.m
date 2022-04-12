%% Load data

cd('~/Desktop/single_cell_manuscript/')
load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed3/reference_bulk_profiles.mat', 'ref', 'r_bulk')
GM12878 = load('data/processed3/GM12878.mat', 'replication_state_filtered', ...
    'percent_replicated_filtered', 'is_included_chr', 'aggregate_S_G1');

GM18507 = load('data/processed3/GM18507.mat', 'replication_state_filtered', ...
    'percent_replicated_filtered', 'is_included_chr', 'aggregate_S_G1');

%% Figure skeleton
figure2 = figure;

set(figure2, 'Position',  [63.5 22.86 18 18.5])

panelA = struct('top', axes('Position', [0.92 16.4438 16.58 1.3078]), ...
'bottom', axes('Position', [0.92 12.564 16.58 3.7055], 'Box', 'off'));

inset = axes('Position', [0.92 11.2126 16.58 1.3514], 'Visible', 'off');

panelA1 = struct('top', axes('Position', [2.1 9.9047 2 1.3078]), ...
    'bottom', axes('Position', [2.1 6.4172 2 3.3131], 'Box', 'off'));

panelA2 = struct('top', axes('Position', [7 9.9047 9.75 1.3078]), ...
    'bottom', axes('Position', [7 6.4172 9.75 3.3131], 'Box', 'off'));

panelB = struct('top', axes('Position', [0.92 3.758 16.58 1.3078]), ...
    'bottom', axes('Position', [0.92 0.75 16.58 2.8336], 'Box', 'off'));

%% Heatmaps

Chr = 2;

params = struct('panel', {panelA panelA1 panelA2 panelB}, 'Chr', 2, ...
    'X', {[100 235] [125 130] [200 230] [100 235]}, ...
    'sample', {GM12878 GM12878 GM12878 GM18507}, ...
    'title', {'GM12878, 10x Genomics' '' '' 'GM18507, DLP+'}, ...
    'inset', {[NaN NaN] [109.5 125.9] [149.3 228.9], [NaN NaN]}, ...
    'name', {'GM12878' 'GM12878' 'GM12878' 'GM18507'});

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
    set(bottom, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [0.5 num_cells+0.5], 'CLim', [2 4], ...
        'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    xlabel(bottom, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    ylabel(bottom, '% Replicated')
    colormap(bottom, [convert_hex(g1_light); convert_hex(s_light{1})])
    yyaxis(bottom, 'right')
    set(bottom, 'YColor', 'k', 'YTick', [])
    ylabel(bottom, [num2str(num_cells) '  cells'])
    
end

xlabel(panelA1.bottom, ['Chr. ' num2str(Chr) ', Mb'])

pos = get(panelA.bottom.XLabel, 'Position');
pos(2) = pos(2) +75;
set(panelA.bottom.XLabel, 'BackgroundColor', 'white', 'Position', pos)

%% Legend for panel A

plot(0, 0, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Bulk-seq', ...
    'Parent', panelA.top)
plot(0, 0, 'Color', s_light{1}, 'LineWidth', 2, 'LineStyle', '-', ...
    'DisplayName', 'Single-Cell S/G1 Aggregate', 'Parent', panelA.top)
legendA = legend(panelA.top);
legendA.ItemTokenSize(1) = 15;
pos = get(panelA.top, 'Position');
set(legendA, 'Units', 'centimeters', 'Orientation', 'horizontal', ...
    'Position', [pos(1) pos(2)+pos(4)+0.1 5.8385 0.4233])

colorbarA = colorbar(panelA.bottom);
set(colorbarA, 'Orientation', 'horizontal', 'Ticks', [],...
    'YAxisLocation', 'top', 'Units', 'centimeters', ...
    'Position', [pos(1)+pos(3)-2.6988 pos(2)+pos(4)+0.1 2.6988 0.4233])
t = ylabel(colorbarA, 'Copy Number', 'FontSize', 7, 'Position', [1.35 0.5], 'VerticalAlignment', 'middle');

cb_labels = axes('Position', colorbarA.Position, 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text([0.25 0.75], [0.5 0.5], {'2' '4'}, 'FontName', 'Arial', 'FontSize', 7, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Parent', cb_labels)

%% Connect insets

set(inset, 'XLim', params(1).X, 'YLim', [0 1], 'Visible', 'off')
uistack(inset, 'bottom')

for p = 2:3
    plot(params(p).X(1) * ones(1, 2), [-2 2], 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off', ...
        'Parent', panelA.top)
    plot(params(p).X(2) * ones(1, 2), [-2 2], 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off', ...
        'Parent', panelA.top)

    plot(params(p).X(1) * ones(1, 2), [0.5 2000], 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off', ...
        'Parent', panelA.bottom)
    plot(params(p).X(2) * ones(1, 2), [0.5 2000], 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off', ...
        'Parent', panelA.bottom)
    
    plot([params(p).X(1) params(p).inset(1)], [1 0], 'k:', 'LineWidth', 0.6, 'Parent', inset);
    plot([params(p).X(2) params(p).inset(2)], [1 0], 'k:', 'LineWidth', 0.6, 'Parent', inset);
end

%% Annotate panels
 
params = struct('panel', {panelA.top panelB.top}, 'text', {'a', 'b'});

for p = 1:2
    text(-0.0397, 1.0638, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/v4/Figure2.pdf')
close
clearvars
