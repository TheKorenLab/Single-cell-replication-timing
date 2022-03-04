load_single_cell_project

load('data/processed/simulated_GM12878.mat',  'true_20kb', 'inferred_state', ...
    'percent_replicated', 'mean_read_depth', 'num_cells', 'p_missed_by_size')
load('data/processed/GM12878.mat', 'single_cell_IRs')
load('data/hg37_genome_metadata.mat', 'genome_windows')

flat_truth = cell2mat(true_20kb);
flat_inference = cell2mat(inferred_state);

is_included = ismember(flat_truth, [2 4]) & ismember(flat_inference, [2 4]);
discordant_windows = sum(flat_truth ~= flat_inference & is_included, 1) ./ sum(is_included, 1);

discordant_IRs = cell(22, 1);
for Chr = 1:22
    
    discordant_IRs{Chr} = false(size(single_cell_IRs{Chr}, 1), num_cells);
    for o = 1:size(single_cell_IRs{Chr}, 1)
        IR_center = find(genome_windows{Chr}(:, 1) < single_cell_IRs{Chr}(o, 4) & ...
            genome_windows{Chr}(:, 2) >= single_cell_IRs{Chr}(o, 4));
        
        missed = true_20kb{Chr}(IR_center, :) == 4 & ...
            inferred_state{Chr}(IR_center, :) == 2;
        extra = true_20kb{Chr}(IR_center, :) == 2 & ...
            inferred_state{Chr}(IR_center, :) == 4;
        
        discordant_IRs{Chr}(o, :) = missed | extra;
    end
end

flat_IRs = cell2mat(discordant_IRs);
discordant_irs_per_cell = sum(flat_IRs, 1);
discordant_cells_per_IR = sum(flat_IRs, 2) ./ num_cells;

IR_rt = cell2mat(single_cell_IRs);
IR_rt = IR_rt(:, 11);


figureS3 = figure;
set(figureS3, 'Units', 'inches', 'Position', [25 9 6.5 4.2])

panelA = axes('Units', 'inches', 'Position', [0.45 2.61 2.9 1.2], 'Box', 'off');
panelB = axes('Units', 'inches', 'Position', [3.5 2.61 2.9 1.2], 'Box', 'off');

w = 1.22;
x = linspace(0.28, 6.4-w, 4);
panelC = axes('Units', 'inches', 'Position', [x(1) 0.5 w 1.2]);
panelD = axes('Units', 'inches', 'Position', [x(2) 0.5 w 1.2]);
panelE = axes('Units', 'inches', 'Position', [x(3) 0.5 w 1.2]);
panelF = axes('Units', 'inches', 'Position', [x(4) 0.5 w 1.2]);

[~, sort_order] = sort(percent_replicated, 'ascend');
[Yticks, YLabels] = get_heatmap_yticks(percent_replicated(sort_order));
 
Chr = 2;
X = [0 76];

params = struct('panel', {panelA panelB}, 'data', {true_20kb inferred_state}, ...
    'title', {'True States' 'Inferred States'});

for p = 1:2
    
    parent = params(p).panel;
    
    imagesc(genome_windows{Chr}(: , 3) ./1e6, 1:num_cells, params(p).data{Chr}(:, sort_order)', ...
        'AlphaData', ~isnan(params(p).data{Chr}(:, sort_order)'), 'Parent', parent)
    set(parent, 'XLim', X, 'YLim', [0.5 num_cells+0.5], 'YDir', 'reverse')
    xlabel(parent, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    title(parent, params(p).title)
    
    if p == 1
        set(parent, 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end))
        ylabel(parent, '% replicated')
    else
        set(parent, 'YTick', [])
    end
    
    colormap(parent, [convert_hex(g1_light); convert_hex(s_light{1})])
    
end

plot(percent_replicated .* 100, discordant_windows .* 100, 'k.', 'MarkerSize', 3, 'Parent', panelC)
set(panelC, 'YLim', [0 5.5])
xlabel(panelC, '% replicated')
ylabel(panelC, '% windows')
title(panelC, 'Effect of S phase')

plot(mean_read_depth, discordant_windows .* 100, 'k.', 'MarkerSize', 3, 'Parent', panelD)
set(panelD, 'YLim', [0 5.5])
xlabel(panelD, 'Reads per 20kb')
ylabel(panelD, '% windows')
title(panelD, 'Effect of coverage')

plot(IR_rt, discordant_cells_per_IR .* 100, 'k.', 'MarkerSize', 3, 'Parent', panelE)
set(panelE, 'XLim', [-1.7 2.15], 'YLim', [0 9])
xlabel(panelE, 'Replication timing')
ylabel(panelE, '% cells')
title(panelE, 'Effect of RT')

x = 1:size(p_missed_by_size, 1);
x = x .* 20e3 ./ 1e6;

plot(x(1:50), p_missed_by_size(1:50), 'k.', 'Parent', panelF)
xlabel(panelF, 'Segment length, Mb')
ylabel(panelF, '% segments missed')
title(panelF, 'Effect of replicon size')
set(panelF, 'YLim', [0 8])

%% Annotate panels

params = struct('panel', {panelA panelB panelC panelD panelE panelF}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f'}, 'x', -0.1007, 'y', 1.1272);

params(2).x = -0.005;

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS3.pdf')
close
clearvars
