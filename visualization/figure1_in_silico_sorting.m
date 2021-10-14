%% Load data

load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed/GM12878.mat', 'subS_fractions', 'aggregate_S_G1')
load('data/processed/reference_bulk_profiles.mat', 'ref', 'r_bulk')

% Scaled MAPD vs. coverage (G1 and S fractions)

fractions = {'G1', 'S'};
mean_coverage = cell(1, 2);
scaled_mapd = cell(1, 2);
is_G1_all = cell(1, 2);
labels = cell(2, 1);
for f = 1:2
    load(['data/intermediate/GM12878-' fractions{f} '-library1-replicate1.mat'], ...
        'scaled_mapd_1Mb', 'mean_coverage_1Mb', 'is_G1')
    
    mean_coverage{f} = mean_coverage_1Mb;
    scaled_mapd{f} = scaled_mapd_1Mb;
    labels{f} = repmat({fractions{f}}, [length(mean_coverage_1Mb) 1]);
    is_G1_all{f} = is_G1;
end

mean_coverage = cell2mat(mean_coverage);
scaled_mapd = cell2mat(scaled_mapd);
is_G1 = cell2mat(is_G1_all);
labels = [labels{1}; labels{2}];

% Raw counts

load('data/hg37_genome_metadata.mat', 'is_mappable')

input_files = {'G1', 'S'};

cells_shown = [1 601];

example_cells.raw = cell(22, 1);
for Chr = 1:22
    coverage_20kb = NaN(size(is_mappable{Chr}, 1), 2);
    for f = 1:length(input_files)
        
        next_library = double(h5read(['data/raw/GM12878-' input_files{f} ...
            '-library1-replicate1.h5'], ['/raw_counts/chr' num2str(Chr)]));
        coverage_20kb(:, f) = next_library(:, cells_shown(f));
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
        example_cells.state{Chr}(:, f) = replication_state{Chr}(:, cells_shown(f));
    end
end

clearvars cells_shown Chr coverage_20kb f fractions input_files is_mappable is_mappable_200kb ...
    next_library replication_state scaled_mapd_1Mb

%% Figure skeleton

figure1 = figure;
set(figure1, 'Position', [25 9 6.5 5.56])

y = [4.91 2.68];

panelA = struct('top', axes('Units', 'inches', 'Position', [0.15 y(1) 1.5 0.3]),...
'bottom', axes('Units', 'inches', 'Position', [0.15 y(1)-0.58 1.5 0.3]));

panelB = struct('left', axes('Units', 'inches', 'Position', [2.3323 y(1)-0.58 1 0.89]), ...
'right', axes('Units', 'inches', 'Position', [3.4167 y(1)-0.58 1 0.89]));

panelC = struct('top', axes('Units', 'inches', 'Position', [4.83 y(1) 1.5 0.3]), ...
'bottom', axes('Units', 'inches', 'Position', [4.83 y(1)-0.58 1.5 0.3]));

panelD = struct('top', axes('Units', 'inches', 'Position', [0.77 y(2) 5.48 0.65]), ...
    'middle', axes('Units', 'inches', 'Position', [0.77 y(2)-0.83 5.48 0.78]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.77 y(2)-2.24 5.48 1.36]));

%% Panel A: GM12878 example cells, raw counts
params = struct('title', {'G1', 'S'}, ...
    'position', {'top', 'bottom'}, 'color', {g1_dark s_dark});

cells_shown = [1 4630];
Chr = 2;
for barcode = 1:2
    
    ymax = prctile(example_cells.raw{Chr}(:, barcode), 99.95) + 10;
    title_string = ['Cell ' num2str(cells_shown(barcode)) ', ' ...
        params(barcode).title];
    
    
    parent = panelA.(params(barcode).position);
    
    plot(genome_windows{Chr}(1:10:end, 3) ./1e6, example_cells.raw{Chr}(:, barcode), ...
        '.', 'Color', params(barcode).color, 'MarkerSize', 2, 'Parent', parent)
    set(parent, 'XLim', [0 76], 'XTick', [], ...
        'YLim', [0 ymax], 'YTick', [])
    title(parent, title_string)
    
    if barcode == 2
        set(parent, 'XTick', 0:25:76)
        xlabel(parent, ['Chr. ' num2str(Chr) ' Coordinate, Mb'])
    end
end
clearvars barcode Chr params ymax title_string parent

%% Panel B: Flow sorting vs. in silico sorting
params = struct('title', {'By FACS', 'In silico'}, ...
    'position', {'left', 'right'}, 'color', {g1_dark s_dark}, 'labels', {'G1', 'S'});
params(1).alpha = 1;
params(2).alpha = 0.5;

for v = 1:2
    
    parent = panelB.(params(v).position);
    
    switch params(v).title
    case 'By FACS'
        sorting = cell(2, 1);
        for fraction = 1:2            
            index = strcmp(labels, params(fraction).labels);
            sorting{fraction} = [mean_coverage(index); scaled_mapd(index)]';
        end
    case 'In silico'
        index = is_G1;
        sorting{1} = [mean_coverage(index); scaled_mapd(index)]';
        index = ~is_G1;
        sorting{2} = [mean_coverage(index); scaled_mapd(index)]';
    end

    for fraction = 1:2
        scatter(sorting{fraction}(:, 1), sorting{fraction}(:, 2), '.', ...
            'MarkerFaceColor', params(fraction).color, ...
            'MarkerFaceAlpha', params(fraction).alpha, ...
            'MarkerEdgeColor', params(fraction).color, ...
            'MarkerEdgeAlpha', params(fraction).alpha, ...
            'HandleVisibility', 'off', 'Parent', parent)
    end
    
    set(parent, 'XLim', [20 270], 'XTick', 50:100:250, 'YLim', [0.8 2.7])
    xlabel(parent, 'Reads per Mb')
    
        
        legend_markers = cell(2, 1);
        for fraction = 1:2
            legend_markers{fraction} = scatter(0, 0, 20, 'filled', 'Parent', parent);
        end
        set(legend_markers{1}, 'MarkerFaceColor', g1_dark, 'DisplayName', 'G1')
        set(legend_markers{2}, 'MarkerFaceColor', s_dark, 'DisplayName', 'S')
        
        if v == 1
            
            set(parent, 'YTick', [1 2])
            ylabel(parent, 'Scaled MAPD')
            
        else
            set(parent, 'YTick', [])
            
        end
      
    title(parent, params(v).title)
    
    panelB.legend.(params(v).position) = legend(parent, 'show');
    x = get(parent, 'Position');
    x = x(1) + 0.01;
    set(panelB.legend.(params(v).position), 'Units', 'inches', 'Orientation', 'vertical', ...
        'Position', [x y(1)-0.05 0.3958 0.3542], 'FontSize', 9)
    panelB.legend.(params(v).position).ItemTokenSize(1) = 5;
end
clearvars fraction index legend_markers params parent sorting v

%% Panel C: GM12878 example cells, replication state
params = struct('title', {'G1', 'S'}, ...
    'position', {'top', 'bottom'}, 'color', {g1_dark s_dark});
Chr = 2;
for barcode = 1:2
    
    title_string = ['Cell ' num2str(cells_shown(barcode)) ', ' ...
        params(barcode).title];
    
    parent = panelC.(params(barcode).position);
    
    plot(genome_windows{Chr}(:, 3) ./1e6, example_cells.state{Chr}(:, barcode), ...
        '.', 'Color', params(barcode).color, 'MarkerSize', 2, 'Parent', parent)
    set(parent, 'XLim', [0 76], 'XTick', [], ...
        'YLim', [1 5], 'YTick', [2 4], 'YTickLabel', {'2N', '4N'})
    title(parent, title_string)
    
    if barcode == 2
        set(parent, 'XTick', 0:25:76)
        xlabel(parent, ['Chr. ' num2str(Chr) ' Coordinate, Mb'])
    end
end
clearvars barcode Chr params ymax title_string parent

%% Panel D: Fractionation

params = struct('title', {'1 fraction', '10 fractions', '100 fractions'}, ...
    'position', {'top', 'middle', 'bottom'});
Chr = 2;
X = [0 76];
cmap = interp1([2 4], [1 1 1; 0 0 0], linspace(2, 4, 64), 'pchip');

set(panelD.top, 'XLim', X, 'XTick', [], 'YTick', [])
plot(ref.GM12878{Chr}(:, 1) ./ 1e6, ref.GM12878{Chr}(:, 2), 'k.', 'HandleVisibility', 'off', ...
    'Parent', panelD.top)
plot(aggregate_S_G1{Chr}(:, 1) ./ 1e6, aggregate_S_G1{Chr}(:, 2), '.', 'Color', s_light{1}, ...
    'HandleVisibility', 'off', 'Parent', panelD.top)
plot(0, 0, 'k', 'LineWidth', 2, 'DisplayName', 'Bulk-seq', 'Parent', panelD.top)
plot(0, 0, 'Color', s_light{1}, 'LineWidth', 2, 'DisplayName', 'Single-Cell S/G1 Aggregate', ...
    'Parent', panelD.top)
ylabel(panelD.top, params(1).title)

yyaxis(panelD.top, 'right')
set(panelD.top, 'YColor', 'k', 'YTick', [])
ylabel(panelD.top, ['r = ' num2str(r_bulk.GM12878, '%0.2f')])

legendD = legend(panelD.top);
legendD.ItemTokenSize(1) = 15;
legend_y = panelD.top.Position(2) + panelD.top.Position(4) + 0.025;
set(legendD, 'Units', 'inches', 'Orientation', 'horizontal', ...
    'Position', [0.77 legend_y 2.7569 0.2014], 'FontSize', 9)

for f = 1:2
    
    parent = panelD.(params(f+1).position);
    
    n_fractions = size(subS_fractions{Chr, f}, 2);
    imagesc(genome_windows{Chr}(:, 3)./1e6, 1:n_fractions, ...
        subS_fractions{Chr, f}', 'AlphaData', ~isnan(subS_fractions{Chr, f}'), 'Parent', parent)
    set(parent, 'YDir', 'reverse', 'YLim', [1 n_fractions], 'Box', 'off')
    
    if f == 1
        set(parent, 'XTick', [])
    else
        xlabel(parent, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    end
    
    set(parent, 'XLim', X, 'YTick', [], 'CLim', [2 4])
    colormap(parent, cmap)
    ylabel(parent, params(f).title)
    
end

colorbarD = colorbar(panelD.bottom);
colorbar_y = panelD.middle.Position(2) + panelD.middle.Position(4) + ...
    panelD.bottom.Position(2) + panelD.bottom.Position(4) - 1.2;
colorbar_y = colorbar_y/2 + panelD.bottom.Position(2);
set(colorbarD, 'Units', 'inches', 'Position', [0.15 colorbar_y-1.2 0.13 1.2], 'Ticks', 2:4)


%% Annotate panels

params = struct('panel', {panelA.top panelB.left panelC.top panelD.top}, ...
    'text', {'a', 'b', 'c', 'd'}, 'x', -0.05, ...
    'y', {1.6 1.1531 1.5 1.3242});
params(2).x = -0.1889;

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/Figure1.pdf')
close
