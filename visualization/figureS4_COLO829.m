load_single_cell_project

load('data/processed/COLO829.mat', 'replication_state_filtered', 'percent_replicated_filtered', ...
    'aggregate_S_G1', 'is_included_chr')
load('data/processed/reference_bulk_profiles.mat', 'ref', 'r_bulk')
load('data/intermediate/COLO829-G1-library1-replicate1.mat', 'scaled_mapd_1Mb', ...
    'mean_coverage_1Mb', 'is_G1')
load('data/hg37_genome_metadata.mat', 'genome_windows')

params = struct('color', {g1_dark s_dark}, 'labels', {'G1', 'S'}, 'alpha', {1 0.5});

%% Figure skeleton
figureS4 = figure;
set(figureS4, 'Units', 'inches', 'Position', [25 9 6.5 3.75])

panelA = axes('Units', 'inches', 'Position', [0.5 2.305 1.3 1.2]);
panelB = axes('Units', 'inches', 'Position', [2.39 2.305 3.83 1.2]);
panelC = axes('Units', 'inches', 'Position', [0.5 0.3889 5.72 1.35]);

%% Panel A

fraction = [is_G1; ~is_G1]';
    
for f = 1:2
    scatter(mean_coverage_1Mb(fraction(:, f)), scaled_mapd_1Mb(fraction(:, f)), '.', ...
        'MarkerFaceColor', params(f).color, 'MarkerFaceAlpha', params(f).alpha, ...
        'MarkerEdgeColor', params(f).color, 'MarkerEdgeAlpha', params(f).alpha, ...
        'HandleVisibility', 'off', 'Parent', panelA)
end
    
Xmax = prctile(mean_coverage_1Mb, 99);
Ymax = prctile(scaled_mapd_1Mb, 99);
set(panelA, 'XLim', [20 Xmax], 'YLim', [0.8 Ymax], 'YTick', [0 1 2])

xlabel(panelA, 'Reads per Mb')
ylabel(panelA, 'Scaled MAPD')
title(panelA, 'COLO-829')

legend_markers = cell(2, 1);
for f = 1:2
    legend_markers{f} = scatter(0, 0, 20, 'filled', 'Parent', panelA);
end
set(legend_markers{1}, 'MarkerFaceColor', g1_dark, 'DisplayName', 'G1')
set(legend_markers{2}, 'MarkerFaceColor', s_dark, 'DisplayName', 'S')

legendA = legend(panelA);
legendA.ItemTokenSize(1) = 5;
pos = get(panelA, 'Position');
set(legendA, 'Units', 'inches', 'Position', [pos(1)+0.05 pos(2)+pos(4)-0.4 0.4 0.35])

%% Panel B

Chr = 2;
X = [0 76];

set(panelB, 'XLim', X, 'YLim', [-2.5 4], 'YTick', [-2 0 2])
plot(0, 0, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Bulk-seq', 'Parent', panelB)
plot(0, 0, 'Color', s_light{1}, 'LineWidth', 2, 'DisplayName', 'Single-Cell S/G1 Aggregate', ...
    'Parent', panelB)
plot(ref.COLO829{Chr}(:, 1) ./ 1e6, ref.COLO829{Chr}(:, 2), '.', 'Color', 'k', ...
    'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', panelB)
plot(aggregate_S_G1{Chr}(:, 1) ./ 1e6, aggregate_S_G1{Chr}(:, 2), '.', ...
    'Color', s_light{1}, 'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', panelB)
xlabel(panelB, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
ylabel(panelB, 'RT')
title(panelB, 'S-phase Cells in the G1 Fraction')
 
legendB = legend(panelB);
legendB.ItemTokenSize(1) = 15;

pos = get(panelB, 'Position');
set(legendB, 'Orientation', 'horizontal', 'Units', 'inches', ...
    'Position', [pos(1)+0.98 pos(2)+pos(4)-0.15 0.43 0.05])

yyaxis(panelB, 'right')
set(panelB, 'YColor', 'k', 'YTick', [])
ylabel(panelB, ['r = ' num2str(r_bulk.COLO829, '%0.2f')])

%% Panel C

num_cells = sum(is_included_chr(Chr, :));
index = percent_replicated_filtered(is_included_chr(Chr, :));
[Yticks, YLabels] = get_heatmap_yticks(index);

r = replication_state_filtered{Chr}(:, is_included_chr(Chr, :));
imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
    'Parent', panelC);
set(panelC, 'YDir', 'reverse', 'XLim', X, 'YLim', [0.5 num_cells+0.5], 'CLim', [2 4], ...
    'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
xlabel(panelC, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'], 'BackgroundColor', 'white')
ylabel(panelC, '% Replicated')
colormap(panelC, [convert_hex(g1_light); convert_hex(s_light{1})])

yyaxis(panelC, 'right')
set(panelC, 'YColor', 'k', 'YTick', [])
ylabel(panelC, [num2str(num_cells) '  cells'])

%% Annotate panels

params = struct('panel', {panelA panelB panelC}, ...
    'text', {'a', 'b', 'c'}, 'x', -0.0634, 'y', 1.0997);

params(1).x = -0.268;
params(2).y = 1.0766;

for p = 1:3
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS4.pdf')
close
clearvars
