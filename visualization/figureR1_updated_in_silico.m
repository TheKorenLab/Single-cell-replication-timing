load_single_cell_project

filename = 'GM12878-G1-library1-replicate1.mat';
load(['data/intermediate/' filename], 'scaled_mapd_1Mb', 'mean_coverage_1Mb');

v1 = load(['data/intermediate/' filename], 'is_G1');
v2 = load(['data/intermediate2/' filename], 'is_G1');

is_G1 = struct('v1', v1.is_G1, 'v2', v2.is_G1);
clearvars v1 v2

coefficients = cell(2, 1);
for v = 1:2
    coefficients{v} = [polyfit(mean_coverage_1Mb(is_G1.(['v' num2str(v)])), ...
        scaled_mapd_1Mb(is_G1.(['v' num2str(v)])), 1);
        polyfit(mean_coverage_1Mb(~is_G1.(['v' num2str(v)])), ...
        scaled_mapd_1Mb(~is_G1.(['v' num2str(v)])), 1)];
end

figureR1 = figure;
set(figureR1, 'Units', 'inches', 'Position', [25 12 3 1.5])

pos = repmat([0.7084 0.776 2.68 2.7], [2 1]);
pos(2, 1) = 4.6;

for v = 1:2
    
    parent = axes('Position', pos(v, :));
    scatter(mean_coverage_1Mb(is_G1.(['v' num2str(v)])), ...
        scaled_mapd_1Mb(is_G1.(['v' num2str(v)])), '.', ...
        'MarkerFaceColor', g1_dark, 'MarkerEdgeColor', g1_dark, ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'Parent', parent, 'HandleVisibility', 'off')
    
    scatter(mean_coverage_1Mb(~is_G1.(['v' num2str(v)])), ...
        scaled_mapd_1Mb(~is_G1.(['v' num2str(v)])), '.', ...
        'MarkerFaceColor', s_dark, 'MarkerEdgeColor', s_dark, ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'Parent', parent, ...
        'HandleVisibility', 'off')
    
    for f = 1:2
        plot(mean_coverage_1Mb, polyval(coefficients{v}(f, :), mean_coverage_1Mb), 'k', ...
            'Parent', parent, 'HandleVisibility', 'off')
    end
    
    title(parent, ['v' num2str(v)])
    xlabel(parent, 'Reads per Mb')
    ylabel(parent, 'Scaled MAPD')
    set(parent, 'XLim', [20 270], 'XTick', 50:100:250, 'YLim', [0.8 2.7], 'YTick', [1 2])

    legend_markers = cell(2, 1);
    for fraction = 1:2
        legend_markers{fraction} = scatter(0, 0, 20, 'filled', 'Parent', parent);
    end
    set(legend_markers{1}, 'MarkerFaceColor', g1_dark, 'DisplayName', 'G1')
    set(legend_markers{2}, 'MarkerFaceColor', s_dark, 'DisplayName', 'S')
    lh = legend(parent, 'Location', 'northwest');
    lh.ItemTokenSize(1) = 5;
    
end

printFigure('out/FigureR1.pdf')
close
clearvars
