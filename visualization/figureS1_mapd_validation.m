load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')

% Scaled MAPD vs. coverage (G1 and S fractions)

fractions = {'G1', 'earlyS', 'S', 'lateS', 'G2'};
mean_coverage = cell(1, 5);
scaled_mapd = cell(1, 5);
is_G1_all = cell(1, 5);
labels = cell(5, 1);
percent_S = NaN(5, 1);
for f = 1:5
    load(['data/intermediate/GM12878-' fractions{f} '-library1-replicate1.mat'], ...
        'scaled_mapd_1Mb', 'mean_coverage_1Mb', 'is_G1')
    
    mean_coverage{f} = mean_coverage_1Mb;
    scaled_mapd{f} = scaled_mapd_1Mb;
    labels{f} = repmat({fractions{f}}, [length(mean_coverage_1Mb) 1]);
    is_G1_all{f} = is_G1;
    percent_S(f) = sum(~is_G1) / length(is_G1);
end

mean_coverage = cell2mat(mean_coverage);
scaled_mapd = cell2mat(scaled_mapd);
is_G1 = cell2mat(is_G1_all);
labels = [labels{1}; labels{2}; labels{3}; labels{4}; labels{5}];
percent_S = percent_S .* 100;

% Raw counts

input_files = {'G1', 'S'};

cells_shown = [963 756; 4 1190];

example_cells.raw = cell(22, 1);
for Chr = 1:22
    coverage_20kb = NaN(size(is_mappable{Chr}, 1), 4);
    for f = 1:length(input_files)
        
        next_library = double(h5read(['data/raw/GM12878-' input_files{f} ...
            '-library1-replicate1.h5'], ['/raw_counts/chr' num2str(Chr)]));
        if f == 1
            index = 1;
        else
            index = 3;
        end
        
        coverage_20kb(:, index:index+1) = next_library(:, cells_shown(f, :));
    end
    coverage_20kb(~is_mappable{Chr}, :) = NaN;
    example_cells.raw{Chr} = aggregate_counts(coverage_20kb, 10);
    is_mappable_200kb = aggregate_counts(is_mappable{Chr}, 10) == 10;
    example_cells.raw{Chr}(~is_mappable_200kb, :) = NaN;

end

% Replication states

example_cells.state = cell(22, 1);
for f = 1:length(input_files)
    load(['data/intermediate/GM12878-' input_files{f} '-library1-replicate1.mat'], ...
        'replication_state')
    for Chr = 1:22
        if f == 1
            index = 1;
        else
            index = 3;
        end
        
        example_cells.state{Chr}(:, index:index+1) = replication_state{Chr}(:, cells_shown(f, :));
    end
end

clearvars Chr coverage_20kb f index input_files is_G1_all is_mappable_200kb ...
    mean_coverage_1Mb next_library replication_state scaled_mapd_1Mb

load('data/processed/GM12878.mat', 'aggregate_G1_fraction')
load('data/processed/reference_bulk_profiles.mat', 'ref')

profiles = cell(22, 1);
for Chr = 1:22
    profiles{Chr}(:, 1) = interp1(ref.GM12878{Chr}(:, 1), ref.GM12878{Chr}(:, 2), ...
        genome_windows{Chr}(:, 3));
    profiles{Chr}(:, 2) = interp1(aggregate_G1_fraction{Chr}(:, 1), ...
        aggregate_G1_fraction{Chr}(:, 2), genome_windows{Chr}(:, 3));
    profiles{Chr}(~is_mappable{Chr}, :) = NaN;
end

profiles = cell2mat(profiles);
r_bulk = corr(profiles(:, 1), profiles(:, 2), 'rows', 'pairwise');

%% Figure
fraction_colors = {'#005A32' '#78C679' '#41AB5D' '#238443' '#005A32'};
cell_colors = {'#E41A1C', '#377EB8', '#984EA3', '#FF7F00'};

fraction_titles = fractions;
fraction_titles{2} = 'Early S';
fraction_titles{4} = 'Late S';

barcodes_shown = [963 756 2328 3514];
fractions_highlighted = {'G1' 'S'};

figureS1 = figure;
set(figureS1, 'Position', [25 9 6.5 7])

panelA = struct('G1', axes('Units', 'inches', 'Position', [0.5 5.684 1 1]), ...
    'earlyS', axes('Units', 'inches', 'Position', [1.6875 5.684 1 1]), ...
    'S', axes('Units', 'inches', 'Position', [2.875 5.684 1 1]), ...
    'lateS', axes('Units', 'inches', 'Position', [4.0625 5.684 1 1]), ...
    'G2', axes('Units', 'inches', 'Position', [5.25 5.684 1 1]));

panelB = struct('left', axes('Units', 'inches', 'Position', [0.5 3.814 1 1]), ...
    'middle_top', axes('Units', 'inches', 'Position', [2 4.249 1.8 0.45]), ...
    'middle_bottom', axes('Units', 'inches', 'Position', [2 3.929 1.8 0.27]), ...
    'right_top', axes('Units', 'inches', 'Position', [4.45 4.249 1.8 0.45]), ...
    'right_bottom', axes('Units', 'inches', 'Position', [4.45 3.929 1.8 0.27]));

panelC = struct('left', axes('Units', 'inches', 'Position', [0.5 1.944 1 1]), ...
    'middle_top', axes('Units', 'inches', 'Position', [2 2.379 1.8 0.45]), ...
    'middle_bottom', axes('Units', 'inches', 'Position', [2 2.059 1.8 0.27]), ...
    'right_top', axes('Units', 'inches', 'Position', [4.45 2.379 1.8 0.45]), ...
    'right_bottom', axes('Units', 'inches', 'Position', [4.45 2.059 1.8 0.27]));

panelD = axes('Units', 'inches', 'Position', [0.5 0.4283 5.75 0.8]);

%% Read depth vs. MAPD (Panel A, Panel B/C left)

panels = {panelA.G1 panelA.earlyS panelA.S panelA.lateS panelA.G2 panelB.left panelC.left};
f = [1:5 1 3];

for p = 1:7
    
    parent = panels{p};
    
    set(parent, 'XLim', [20 450], 'XTick', 50:150:450, 'YLim', [0.8 3.8], 'YTick', 1:3)
    xlabel(parent, 'Reads per Mb')
    title(parent, fraction_titles{f(p)})
    
    if ismember(p, [1 6 7])
        ylabel(parent, 'Scaled MAPD')
    else
        set(parent, 'YTick', [])
    end
    
    scatter(mean_coverage, scaled_mapd, '.', ...
        'MarkerFaceColor', convert_hex('#bdbdbd'), 'MarkerEdgeColor', convert_hex('#bdbdbd'), ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Parent', parent)
    
    is_in_fraction = strcmp(labels, fractions{f(p)});
    
    if p <= 5
        scatter(mean_coverage(is_in_fraction), scaled_mapd(is_in_fraction), '.', ...
            'MarkerFaceColor', convert_hex(fraction_colors{p}), ...
            'MarkerEdgeColor', convert_hex(fraction_colors{p}), ...
            'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Parent', parent)
        text(208, 3.285, [num2str(percent_S(p), '%0.1f') '% of cells'  ...
            newline 'replicating'], 'FontSize', 9, 'HorizontalAlignment', 'center', 'Parent', parent)
    else
        is_highlighted = barcodes_shown(is_in_fraction(barcodes_shown));
        for b = 1:length(is_highlighted)
            c = find(barcodes_shown == is_highlighted(b));
            scatter(mean_coverage(is_highlighted(b)), scaled_mapd(is_highlighted(b)), ...
                12, 'MarkerFaceColor', convert_hex(cell_colors{c}), ...
                'MarkerEdgeColor', convert_hex(cell_colors{c}), ...
                'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Parent', parent)
        end
        
    end
end

%% Example cells (Panels B & C, middle & right)

panels = {'middle' 'right' 'middle' 'right'};
Chr = 2;
X = [0 76];

for b = 1:4

    if b <= 2
        p = panelB;
    else
        p = panelC;
    end
    
    top = p.([panels{b} '_top']);
    ymax = prctile(example_cells.raw{Chr}(:, b), 99.95) + 10;
    set(top, 'XTick', [], 'XLim', X, 'YLim', [0 ymax], 'YTick', [])
    plot(genome_windows{Chr}(1:10:end, 3) ./ 1e6, ...
        example_cells.raw{Chr}(:, b), '.', ...
        'MarkerSize', 2, 'Color', cell_colors{b}, 'Parent', top)
    ylabel(top, 'Count')
    title(top, ['Cell ' num2str(barcodes_shown(b)) ', MAPD: ' ...
        num2str(scaled_mapd(barcodes_shown(b)), '%0.2f')])
    
    bottom = p.([panels{b} '_bottom']);
     set(bottom, 'XLim', X, 'YLim', [1 5], 'YTick', [2 4], 'YTickLabel', {'2N', '4N'})
    plot(genome_windows{Chr}(:, 3) ./ 1e6, example_cells.state{Chr}(:, b), '.', ...
        'MarkerSize', 2, 'Color', cell_colors{b}, 'Parent', bottom)
    xlabel(bottom, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    
end

%% Panel D

set(panelD, 'XLim', X, 'YLim', [-2.5 4], 'YTick', [-2 0 2])
plot(0, 0, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Bulk-seq', 'Parent', panelD)
plot(0, 0, 'Color', s_light{1}, 'LineWidth', 2, 'DisplayName', 'Single-Cell S/G1 Aggregate', ...
    'Parent', panelD)
plot(ref.GM12878{Chr}(:, 1) ./ 1e6, ref.GM12878{Chr}(:, 2), '.', 'Color', 'k', ...
    'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', panelD)
plot(aggregate_G1_fraction{Chr}(:, 1) ./ 1e6, aggregate_G1_fraction{Chr}(:, 2), '.', ...
    'Color', s_light{1}, 'MarkerSize', 4, 'HandleVisibility', 'off', 'Parent', panelD)
xlabel(panelD, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
ylabel(panelD, 'Rep. Timing')
title(panelD, 'S-phase Cells in the G1 Fraction')
 
legendD = legend(panelD);
set(legendD, 'Orientation', 'horizontal', 'FontSize', 9, 'Position', [0.0823 0.1468 0.4241 0.0288])
legendD.ItemTokenSize(1) = 15;

yyaxis(panelD, 'right')
set(panelD, 'YColor', 'k', 'YTick', [])
ylabel(panelD, ['r = ' num2str(r_bulk, '%0.2f')])

%% Annotate panels

params = struct('panel', {panelA.G1 panelB.left panelC.left panelD}, ...
    'text', {'a', 'b', 'c', 'd'}, 'x', -0.375, 'y', 1.1806);
params(4).x = -0.0628;
params(4).y = 1.2435;

for p = 1:4
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS1.pdf')
close
