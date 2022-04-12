cd('~/Desktop/single_cell_manuscript/')
load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed3/GM12878.mat', 'single_cell_IRs', 'replication_state_filtered', ...
    'percent_replicated_filtered', 'is_included_chr', 'aggregate_S_G1', 'IR_range', ...
    'firing_order_frequency')

for Chr = 1:22
    single_cell_IRs{Chr} = [Chr * ones(size(single_cell_IRs{Chr}, 1), 1) ...
        single_cell_IRs{Chr} IR_range{Chr}];
end

flat_IR_list = cell2mat(single_cell_IRs);
flat_IR_list = sortrows(flat_IR_list, 12, 'descend');
flat_IR_list(:, 14:15) = flat_IR_list(:, 14:15) .* 100;
flat_IR_list(:, 17) = flat_IR_list(:, 15) - flat_IR_list(:, 14);
flat_IR_list = flat_IR_list(:, [1:15 17 16]);
flat_IR_list(isnan(flat_IR_list(:, 14)), 17) = false;

flat_IR_list(:, 18:19) = NaN;
p = percent_replicated_filtered .* 100;
for row = 1:size(flat_IR_list, 1)
    if isnan(flat_IR_list(row, 14))
        continue
    end
    flat_IR_list(row, 18) = find(p == flat_IR_list(row, 14));
    
    if isnan(flat_IR_list(row, 15))
        continue
    end
    flat_IR_list(row, 19) = find(p == flat_IR_list(row, 15));
end
clearvars p row

Xticks = flip(interp1(flat_IR_list(:, 12), 1:size(flat_IR_list, 1), [-1.25 -1:0.5:1.5]));

num_unexpected_origins = sum(firing_order_frequency(:, 3:4), 2);
f = polyfit(flat_IR_list(:, 12), num_unexpected_origins, 2);
trend_line = polyval(f, flat_IR_list(:, 12));
trend_line(trend_line < 0 | trend_line > 100) = NaN;

sc_colormap = [convert_hex(g1_light); convert_hex(s_light{1})];
cmap = interp1([0 100], [215 48 39; 145 191 219]./255, linspace(0, 100, 128), 'pchip');
bar_colors = {'#DE2C26', '#3182BD', '#636362', '#BDBDBD'};

%% Figure skeleton

figure4 = figure;

set(figure4, 'Position', [63.5 22.86 18 18.5])

panelA = struct('top', axes('Position', [6.35 16.9 5.75 1.2]), ...
    'bottom', axes('Position', [6.35 14.75 5.75 2.1], 'Box', 'off'));

panelB = axes('Position', [0.92 15.3 2.8 2]);

panelC = axes('Position', [0.92 9.6 6.25 3.4]);

panelD = axes('Position', [2.045 4.68 4 2.5]);

panelE = axes('Position', [14.95 15.3 2.8 2]);

panelF = axes('Position', [12.35 10.3 5.4 3.1]);

panelG = axes('Position', [12.35 5.53 5.4 3.1]);

panelH = axes('Position', [0.92 0.7 16.86 2.94]);

%% Panel A

params = struct('Chr', 1, 'X', [74 79.5], 'Y', [-1.5 1.7], 'IR_shown', {208 210 215}, ...
    'IR_X', NaN, 'IR_earliest', NaN, 'IR_latest', NaN, 'firing_stats', NaN);

for IR = 1:size(params, 2)
    index = flat_IR_list(:, 1) == params(IR).Chr & ...
        flat_IR_list(:, 2) == params(IR).IR_shown;
    params(IR).IR_X = flat_IR_list(index, 5) ./ 1e6;
    params(IR).IR_earliest = flat_IR_list(index, 18);
    params(IR).IR_latest = flat_IR_list(index, 19);
    
    params(IR).firing_stats = firing_order_frequency(index, :);
end

plot(aggregate_S_G1{params(1).Chr}(:, 1)./1e6, aggregate_S_G1{params(1).Chr}(:, 2), '.', ...
    'Color', s_dark, 'Parent', panelA.top)
set(panelA.top, 'XLim', params(1).X, 'YLim', params(1).Y, 'XTick', [])
ylabel(panelA.top, 'RT')
title(panelA.top, ['Chromosome ' num2str(params(1).Chr)])

num_cells = sum(is_included_chr(params(1).Chr, :));
[Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(params(1).Chr, :)));

r = replication_state_filtered{params(1).Chr}(:, is_included_chr(params(1).Chr, :));
imagesc(genome_windows{params(1).Chr}(:, 3)./1e6, 1:num_cells, ...
    r', 'AlphaData', ~isnan(r)', 'Parent', panelA.bottom)
set(panelA.bottom, 'YDir', 'reverse', 'XLim', params(1).X, 'YLim', [0.5 num_cells+0.5], ...
    'YTick', Yticks(2:2:end), ...
    'YTickLabel', YLabels(2:2:end), 'Box', 'off')
colormap(panelA.bottom, sc_colormap)

xlabel(panelA.bottom, ['Chromosome ' num2str(params(1).Chr) ' Coordinate, Mb'], ...
    'BackgroundColor', 'white')
xl = get(panelA.bottom, 'XLabel');
set(xl, 'Units', 'centimeters');
pos = xl.Position;
pos(2) = -0.5;
set(xl, 'Position', pos)

ylabel(panelA.bottom, '% Replicated')

for IR = 1:size(params, 2)
    plot(params(IR).IR_X * ones(1,2), params(IR).Y, 'k', 'LineWidth', 0.7, 'Parent', panelA.top)
    plot(params(IR).IR_X * ones(1,2), [1 num_cells], 'k', 'LineWidth', 0.7, 'Parent', panelA.bottom)
    
    plot(params(IR).IR_X * ones(1, 2), [params(IR).IR_earliest params(IR).IR_latest], '*', ...
        'Color', '#8856A7', 'Parent', panelA.bottom)
end

yyaxis(panelA.bottom, 'right')
set(panelA.bottom, 'YColor', 'k', 'YTick', [])
ylabel(panelA.bottom, [num2str(num_cells) ' cells'])

text_params = struct('X', {-0.3289 1.1368}, 'Y', 0.8208, 'text', {'Order' 'Timing'});

for p = 1:2
   text(text_params(p).X, text_params(p).Y, text_params(p).text, 'FontSize', 7, ...
        'Parent', panelA.bottom, 'Units', 'normalized')
end

%% Panel B

example_firing_stats = NaN(3, 4);
example_x = NaN(3, 1);
for IR = 1:size(params, 2)
    example_firing_stats(IR, :) = params(IR).firing_stats;
    example_x(IR) = params(IR).IR_X;
end

b = bar(flip(example_firing_stats, 2), 'stacked', 'EdgeColor', 'none', 'Parent', panelB);
for row = 1:4
    b(row).FaceColor = bar_colors{row};
end

set(panelB, 'XLim', [0.4 3.6], 'YLim', [0 100], 'YTick', [0 25 75], ...
    'XTick', 1:3, 'XTickLabel', num2str(example_x, '%.1f'))

xlabel(panelB, ['Chr' num2str(params(1).Chr) ', Mb'], 'BackgroundColor', 'white')
xl = get(panelB, 'XLabel');
set(xl, 'Units', 'centimeters');
pos = xl.Position;
pos(2) = -0.5;
set(xl, 'Position', pos)

ylabel(panelB, '% of cells')
title(panelB, 'Firing Order')

legendB = legend(panelB, b([4 3 2 1]), {'Fired as predicted', 'Unfired as predicted', ...
    'Delayed firing', 'Premature firing'});
legendB.ItemTokenSize(1) = 10;
set(legendB, 'Units', 'centimeters', 'Position', [7.27 11.6771 3 1.3229])

%% Panel C

b = bar(flip(firing_order_frequency, 2), 1, 'stacked', 'EdgeColor', 'none', 'Parent', panelC);
for row = 1:4
    b(row).FaceColor = bar_colors{row};
end
set(panelC, 'XLim', [1 size(firing_order_frequency, 1)], 'YLim', [0 100], 'YTick', [0 25 75], ...
    'XTick', Xticks, 'XTickLabel', flip([-1.25 -1:0.5:1.5]))
xlabel(panelC, 'Aggregate replication timing')
ylabel(panelC, '% of cells')
title(panelC, 'Predictability of IR Firing Order')

%% Panel D

X = prctile(flat_IR_list(:, 12), [0 100]);
Y = [0 max(num_unexpected_origins(~isoutlier(num_unexpected_origins))) + 5];

plot(flat_IR_list(:, 12), num_unexpected_origins, '.', 'MarkerSize', 4, 'Color', '#BDBDBD', ...
    'Parent', panelD)
plot(flat_IR_list(:, 12), trend_line, 'Color', '#8856A7', 'LineWidth', 2, 'Parent', panelD)
set(panelD, 'XLim', X, 'YLim', Y, 'XDir', 'reverse')
xlabel(panelD, 'Aggregate replication timing')
ylabel(panelD, '% of cells')
title(panelD, 'Out-of-Order IR Firing')

text_params = struct('X', 0.77, 'Y', {1.5521 1.4326 1.3132}, ...
    'Color', {bar_colors{2} 'k' bar_colors{1}}, 'text', {'Delayed', '+', 'Premature'});

for p = 1:size(text_params, 2)
    t=text(text_params(p).X, text_params(p).Y, text_params(p).text, 'Color', text_params(p).Color, ...
        'FontSize', 7, 'HorizontalAlignment', 'center', 'Parent', panelD, 'Units', 'normalized');
end

%% Panel E

for IR = 1:size(params, 2)
    plot(IR * ones(1, 2), [params(IR).IR_earliest params(IR).IR_latest], 'Color', '#8856A7', ...
        'LineWidth', 0.7, 'Marker', '*', 'Parent', panelE)
end
set(panelE, 'YDir', 'reverse', 'XLim', [0.4 3.6],  'YLim', [-100 num_cells + 100], ...
    'XTick', 1:3, 'XTickLabel', num2str(example_x, '%.1f'), ...
    'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end))

xlabel(panelE, ['Chr' num2str(params(1).Chr) ', Mb'])
xl = get(panelE, 'XLabel');
set(xl, 'Units', 'centimeters');
pos = xl.Position;
pos(2) = -0.5;
set(xl, 'Position', pos)

ylabel(panelE, '% Replicated')
title(panelE, 'Firing Range')

%% Panels F & G

params = struct('panel', {panelF panelG}, 'col', {14 15}, ...
    'title', {'Earliest Firing Time per IR' 'Latest Firing Time Per IR'});

index = logical(flat_IR_list(:, 17));

for p = 1:2

    parent = params(p).panel;

    plot(flat_IR_list(:, 12), flat_IR_list(:, params(p).col), '.', 'Color', '#BDBDBD', ...
        'MarkerSize', 4, 'Parent', parent)
    plot(flat_IR_list(index, 12), flat_IR_list(index, params(p).col), '.', ...
        'Color', 'k', 'MarkerSize', 4, 'Parent', parent)
    set(parent, 'XDir', 'reverse', 'YDir', 'reverse', 'YLim', [0 100], 'XLim', X, 'YTick', 20:20:80)
    xlabel(parent, 'Aggregate replication timing')
    ylabel(parent, '% Replicated')
    title(parent, params(p).title)
    
end

%% Panel H

c = interp1(linspace(0, 100, 128), cmap, flat_IR_list(:, 16));
for o = 1:size(flat_IR_list, 1)
    if ~index(o)
        continue
    end
    plot(o * ones(1, 2), flat_IR_list(o, 14:15), 'Color', c(o, :), 'LineWidth', 0.6, ...
        'Parent', panelH)
end
set(panelH, 'YDir', 'reverse', 'XLim', [1 size(flat_IR_list, 1)] + [-5 5], 'CLim', [0 100], ...
    'XTick', Xticks, 'XTickLabel', flip([-1.25 -1:0.5:1.5]), 'YTick', 0:20:100)
xlabel(panelH, 'Aggregate replication timing')
ylabel(panelH, '% Replicated')
title(panelH, 'Range of Firing Time per IR')
colormap(panelH, cmap)

colorbarH = colorbar(panelH);
pos = get(panelH, 'Position');
set(colorbarH, 'Units', 'centimeters', 'Orientation', 'horizontal', ...
    'Position', [pos(1)+pos(3)-2.55 pos(2)+pos(4)+0.08 2.55 0.4], 'Ticks', 20:20:90)

ylabel(colorbarH, {'% of' 'S phase'}, 'FontSize', 7, 'Position', [-21.7 0.5], ...
    'VerticalAlignment', 'middle')

%% Arrows

% A->B, A->E
params = struct('X', {0.2959 0.7123}, 'dir', {-1 1});

for p = 1:2
    annotation('arrow', params(p).X + [0 params(p).dir * 0.0748], 0.8772 * ones(1, 2), ...
        'Color', 'k', 'LineWidth', 0.6);
end


% B->C, E->F
params = struct('X', {0.1723 0.8323}, 'Y', {0.7953 0.8170}, 'dir', {1 -1});

for p = 1:2
    annotation('arrow', params(p).X + [0 params(p).dir * 0.032], params(p).Y + [0 -0.0655], ...
        'Color', 'k', 'LineWidth', 0.6);
end

% C->D
annotation('arrow', [0.2264 0.2264], [0.4645 0.4121], 'Color', 'k', 'LineWidth', 0.6)

% F->H
X = [0.5269 0.5911 0.6376];
Y = [0.233 0.3831 0.5033 0.6411];

params = struct('X', {[2 3], [2 3], [2 2], [1 2], [1 1]}, ...
    'Y', {[4 4], [2 2], [2 4], [3 3], [3 1]}, 'style', 'line');
params(5).style = 'arrow';

for p = 1:5
   annotation(params(p).style, X(params(p).X), Y(params(p).Y), 'Color', 'k', 'LineWidth', 0.7)
end

%% Annotate panels

params = struct('panel', {panelA.top panelB panelC panelD panelE panelF panelG panelH}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}, ...
    'x', {-0.0797 -0.2482 -0.1043, -0.0833 -0.2011 -0.0625 -0.0625 -0.0361}, ...
    'y', {1.1182 1.1528 1.1226 1.1806 1.094 1.1272 1.1272 1.1064});

for p = 1:8
   t=text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/v4/Figure4.pdf')
close
clearvars
