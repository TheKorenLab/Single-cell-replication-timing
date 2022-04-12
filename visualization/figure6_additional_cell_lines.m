cd('~/Desktop/single_cell_manuscript/')
load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed3/cell_line_pca.mat', 'pcs', 'pca_labels')
load('data/processed3/reference_bulk_profiles.mat', 'ref', 'r_bulk')

data = struct;
for sample = 1:9
    data.(samples{sample}) = load(['data/processed3/' samples{sample} '.mat'], ...
        'replication_state_filtered', 'percent_replicated_filtered', 'is_included_chr', ...
        'aggregate_S_G1', 'single_cell_IRs', 'late_IR_class', 'replication_tracks');
end

cell_line_colors = {'#41AE76','#238B45', '#005824', '#74A9CF', '#2B8CBE', '#045A8D', '#F768A1', ...
    '#C51B8A', '#7A0177'};

%% Figure skeleton

figure6 = figure;
set(figure6, 'Position', [63.5 22.86 18 18.5])

panelA = struct('top', axes('Position', [0.92 17.1 7.4 1]), ...
    'bottom', axes('Position', [0.92 15 7.4 2]));

panelB = struct('top', axes('Position', [10.1 17.1 7.4 1]), ...
    'bottom', axes('Position', [10.1 15 7.4 2]));

panelC = struct('bulk', axes('Position', [0.92 12.3 4.5 1]), ...
    'sc1', axes(), 'sc2', axes(), 'sc3', axes(), 'sc4', axes(), 'sc5', axes(), 'sc6', axes());

panelD = axes('Position', [11.5 9.8 5.1 3.5]);

panelE = struct('left', axes(), 'middle', axes(), 'middle_inset', axes(), 'right', axes());

insetE = axes('Position', [12.1 5.2 0.65 3.1], 'Visible', 'off');

panelF = struct('left', axes(), 'middle', axes(), 'middle_inset', axes(), 'right', axes());

insetF = axes('Position', [12.1 0.75 0.65 3.1], 'Visible', 'off');

panels = {panelC.sc1 panelC.sc2 panelC.sc3 panelC.sc4 panelC.sc5 panelC.sc6};
y = linspace(10.368, 0.75, 6);

for p = 1:6
    set(panels{p}, 'Units', 'centimeters', 'Position', [0.92 y(p) 4.5 1.85])
end

panels = {panelE.left panelE.middle panelE.middle_inset panelE.right ...
    panelF.left panelF.middle panelF.middle_inset panelF.right};

x = [7.25 10.25 12.75 15.5];
y = [5.2 0.75];

pos = cell(4, 2);
for xi = 1:4
    for yi = 1:2
        pos{xi, yi} = [x(xi) y(yi) 1.85 3.1];
    end
end
pos = cell2mat(pos(:));

for p = 1:8
    set(panels{p}, 'Units', 'centimeters', 'Position', pos(p, :))
    
end

%% Panel A/B

params = struct('panel', {panelA panelB}, 'Chr', {4 5}, 'X', {[75 120] [1 45]}, ...
    'sample', {data.H7 data.MCF7}, 'bulk', {ref.H7 ref.MCF7}, 'r_bulk', {r_bulk.H7 r_bulk.MCF7}, ...
    'color', {s_light{2} s_light{3}}, 'title', {'H7-hESC' 'MCF-7 Breast Cancer'});

for p = 1:2
    
    plot(params(p).bulk{params(p).Chr}(:, 1) ./ 1e6, params(p).bulk{params(p).Chr}(:, 2), 'k.', ...
        'Parent', params(p).panel.top)
    plot(params(p).sample.aggregate_S_G1{params(p).Chr}(:, 1) ./ 1e6, ...
        params(p).sample.aggregate_S_G1{params(p).Chr}(:, 2), '.', 'Color', params(p).color, ...
        'Parent', params(p).panel.top)
    set(params(p).panel.top, 'XLim', params(p).X, 'XTick', [])
    ylabel(params(p).panel.top, 'RT')
    title(params(p).panel.top, params(p).title)
    
    yyaxis(params(p).panel.top, 'right')
    set(params(p).panel.top, 'YColor', 'k', 'YTick', [])
    ylabel(params(p).panel.top, ['r = ' num2str(params(p).r_bulk, '%0.2f')])
    
    index = params(p).sample.is_included_chr(params(p).Chr, :);
    num_cells = sum(index);
    [Yticks, YLabels] = get_heatmap_yticks(params(p).sample.percent_replicated_filtered(index));
    r = params(p).sample.replication_state_filtered{params(p).Chr}(:, index);
    imagesc(genome_windows{params(p).Chr}(:, 3) ./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
        'Parent', params(p).panel.bottom)
    set(params(p).panel.bottom, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [0.5 num_cells+0.5], ...
        'CLim', [2 4], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(params(p).panel.bottom, [convert_hex(g1_light); convert_hex(params(p).color)])
    xlabel(params(p).panel.bottom, ['Chromosome ' num2str(params(p).Chr) ' Coordinate, Mb'])
    ylabel(params(p).panel.bottom, '% Replicated')
    
    yyaxis(params(p).panel.bottom, 'right')
    set(params(p).panel.bottom, 'YColor', 'k', 'YTick', [])
    ylabel(params(p).panel.bottom, [num2str(num_cells) ' cells'])
    
end

%% Panel C

params = struct('Chr', 4, 'X', [58 70]);

for sample = 1:6
    plot(ref.(samples{sample}){params.Chr}(:, 1) ./ 1e6, ref.(samples{sample}){params.Chr}(:, 2), ...
        '.', 'Color', cell_line_colors{sample}, 'Parent', panelC.bulk)
end
set(panelC.bulk, 'XLim', params(1).X, 'XTick', [], 'YLim', [-2.5 2.5], 'YTick', [-1 1])
ylabel(panelC.bulk, 'RT')
title(panelC.bulk, 'Cell-Type Specificity')

params = struct('Chr', 4, 'X', [58 70], ...
    'sample', {data.GM12878 data.GM12891 data.GM12892 data.H1 data.H7 data.H9}, ...
    'color', {s_light{1} s_light{1} s_light{1} s_light{2} s_light{2} s_light{2}}, ...
    'panel', {panelC.sc1 panelC.sc2 panelC.sc3 panelC.sc4 panelC.sc5 panelC.sc6}, ...
    'title', {'GM12878' 'GM12891' 'GM12892' 'H1', 'H7', 'H9'});

for p = 1:6
    
    index = params(p).sample.is_included_chr(params(p).Chr, :);
    num_cells = sum(index);
    [Yticks, YLabels] = get_heatmap_yticks(params(p).sample.percent_replicated_filtered(index));
    r = params(p).sample.replication_state_filtered{params(p).Chr}(:, index);
    imagesc(genome_windows{params(p).Chr}(:, 3) ./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
        'Parent', params(p).panel)
    set(params(p).panel, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [0.5 num_cells+0.5], ...
        'CLim', [2 4], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(params(p).panel, [convert_hex(g1_light); convert_hex(params(p).color)])
    ylabel(params(p).panel, '% Rep.')
    
    if p <= 5
        set(params(p).panel, 'XTick', [])
    else
        xlabel(params(p).panel, ['Chromosome ' num2str(params(p).Chr) ' Coordinate, Mb'])
    end
    
    yyaxis(params(p).panel, 'right')
    set(params(p).panel, 'YColor', 'k', 'YTick', [])
    ylabel(params(p).panel, params(p).title)
end

%% Panel D

for sample = 1:9
    index = strcmp(pca_labels, samples{sample});
    scatter(pcs(index, 1), pcs(index, 2), 5, convert_hex(cell_line_colors{sample}), 'filled', ...
        'HandleVisibility', 'off', 'Parent', panelD)
    plot(-500, 0, '.', 'Color', convert_hex(cell_line_colors{sample}), 'MarkerSize', 15, ...
        'Parent', panelD);
end
set(panelD, 'XLim', [-230 230], 'YLim', [-110 110])
xlabel(panelD, 'Principal Component 1')
ylabel(panelD, 'Principal Component 2')
title(panelD, 'Replication Trajectories')

legendD = legend(panelD, cell_line_names(1:9));
set(legendD, 'Units', 'inches', 'Position', [2.85 3.9 0.7222 1.1111])
legendD.ItemTokenSize(1) = 10;

%% Panel E

params = struct('panel', {panelE.left panelE.middle panelE.right ...
    panelF.left panelF.middle panelF.right}, 'Chr', {3 2 10 6 5 2}, ...
    'sample', {data.H7 data.H7 data.H7 data.MCF7 data.MCF7 data.MCF7}, ...
    'IR_shown', {312 8 21 58 28 53}, 'inset', [], 'inset_Y', NaN, ...
    'color', {s_light{2} s_light{2} s_light{2} s_light{3} s_light{3} s_light{3}});

params(2).inset = insetE;
params(5).inset = insetF;

for p = 1:size(params, 2)
    
    parent = params(p).panel;
    
    X = params(p).sample.single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4)/1e6 + [-.5 .5];
    
    index = params(p).sample.is_included_chr(params(p).Chr, :);
    num_cells = sum(index);
    [Yticks, YLabels] = get_heatmap_yticks(params(p).sample.percent_replicated_filtered(index));
    
    r = params(p).sample.replication_state_filtered{params(p).Chr}(:, index);
    imagesc(genome_windows{params(p).Chr}(:, 3)./ 1e6, 1:num_cells, r', 'AlphaData', r', ...
        'Parent', parent)
    set(parent, 'YDir', 'reverse', 'XLim', X, 'YLim', [0.5 num_cells+0.5], 'YTick', Yticks(2:2:end), ...
        'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(params(p).panel, [convert_hex(g1_light); convert_hex(params(p).color)])
    plot(params(p).sample.single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 .* ...
        ones(1,2), [0.5 num_cells+0.5], 'k', 'LineWidth', 0.6, 'Parent', parent)
    xlabel(parent, ['Chr' num2str(params(p).Chr) ', Mb'])
    
    if ismember(p, [1 4])
        ylabel(parent, '% Replicated')
    end
    
    if params(p).sample.late_IR_class{params(p).Chr}(params(p).IR_shown, 1)
        title(parent, 'Const. Late');
    elseif params(p).sample.late_IR_class{params(p).Chr}(params(p).IR_shown, 2)
        title(parent, 'Early+Rare');
    elseif params(p).sample.late_IR_class{params(p).Chr}(params(p).IR_shown, 3)
        title(parent, 'Throughout S');
    end
    
    if ~isempty(params(p).inset)
        set(params(p).inset, 'YDir', 'reverse', 'XLim', [0 1], 'YLim', [0.5 num_cells+0.5], ...
            'Visible', 'off')
        params(p).inset_Y = num_cells;
    end
    
end

% Zoom in
params(2).panel = panelE.middle_inset;
params(5).panel = panelF.middle_inset;

for p = [2 5]
    parent = params(p).panel;
    
    X = params(p).sample.single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4)/1e6 + [-0.5 0.5];
    
    index = params(p).sample.replication_tracks{params(p).Chr}(:, 8) == params(p).IR_shown;
    barcodes = params(p).sample.replication_tracks{params(p).Chr}(index, 1);
    
    YLabels = params(p).sample.percent_replicated_filtered(barcodes) .* 100;
    
    r_tmp = params(p).sample.replication_state_filtered{params(p).Chr}(:, params(p).sample.is_included_chr(params(p).Chr, :));
    
    barcodes_shown = find(params(p).sample.is_included_chr(params(p).Chr, :));
    
    num_cells = length(barcodes);
    
    pos = get(parent, 'Position');
    pos(4) = 0.6 * num_cells;
    pos(2) = pos(2) + (0.6 * (5 - num_cells));
    set(parent, 'Position', pos)
    
    col = (num_cells * 11) + (2 * (num_cells-1));
    r = NaN(size(r_tmp, 1), col);
    Y = NaN(num_cells, 1);
    counter = 1;
    for b = 1:num_cells
        barcode = find(barcodes_shown == barcodes(b));
        Y(b) = barcode;
        barcode = barcode +  [-5 5];
        if barcode(1) < 1
            barcode(1) = 1;
            r(:, counter+1:counter+10) = r_tmp(:, barcode(1):barcode(2));
        else
            r(:, counter:counter+10) = r_tmp(:, barcode(1):barcode(2));
        end
        
        if b < size(barcodes, 1)
            r(:, counter+11:counter+12) = NaN;
            counter = counter + 13;
        end
    end
    Y = Y([1 end]);
    
    imagesc(genome_windows{params(p).Chr}(:, 3) ./ 1e6, 1:col, r', 'AlphaData', ~isnan(r'), ...
        'Parent', parent)
    set(parent, 'XLim', X, 'YLim', [1 col], 'YDir', 'reverse', 'YTick', [6 19 32 45 58], ...
        'YTickLabel', round(YLabels, 0), 'Box', 'off')
    colormap(params(p).panel, [convert_hex(g1_light); convert_hex(params(p).color)])
    
    plot(params(p).sample.single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 * ones(1,2), [1 63], ...
        'k', 'LineWidth', 0.6, 'Parent', parent)
    
    Y(:, 2) = [0 params(p).inset_Y / 3.1 * pos(4)];
    if p == 2
        Y(:, 2) = Y(:, 2) + 10;
    else
        Y(:, 2) = Y(:, 2) + 4;
    end
    
    for l = 1:2
      plot([0 1], Y(l, :), 'k:', 'LineWidth', 0.6, 'Parent', params(p).inset)
    end
    
end

%% Annotate panels

params = struct('panel', {panelA.top panelB.top panelC.bulk panelD panelE.left panelF.left}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f'}, ...
    'x', {-0.0857 -0.0857 -0.1392 -0.1 -0.25 -0.25}, ...
    'y', {1.1441 1.1441 1.2 1.1 1.05 1.05});

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/v4/Figure6.pdf')
close
clearvars
