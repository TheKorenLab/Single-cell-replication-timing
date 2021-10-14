load_single_cell_project

load('data/processed/GM12878.mat', 'p_cells_fired', 'single_cell_IRs', ...
    'replication_state_filtered', 'late_IR_class', 'aggregate_S_G1', 'percent_replicated_filtered', ...
    'is_included_chr', 'replication_tracks')
load('data/hg37_genome_metadata.mat', 'genome_windows')

sc_colormap = [convert_hex(g1_light); convert_hex(s_light{1})];

flat_p_fired = cell2mat(p_cells_fired);
flat_IRs = cell2mat(single_cell_IRs);
flat_classes = cell2mat(late_IR_class);

flat_classes = flat_classes(:, [3 2 1]);

%% Figure skeleton

figure5 = figure;
set(figure5, 'Position', [25 9 6.5 8.5])

panelA = struct('top', axes('Units', 'inches', 'Position', [0.45 7.7 4.1 0.5]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.45 6.6 4.1 1.05]));

panelB = struct('top_left', axes('Units', 'inches', 'Position', [0.45 4.46 0.7 1.2]), ...
    'top_right', axes('Units', 'inches', 'Position', [1.59 4.46 0.7 1.2]), ...
    'bottom_left', axes('Units', 'inches', 'Position', [0.45 2.45 0.7 1.2]), ...
    'bottom_right', axes('Units', 'inches', 'Position', [1.59 2.45 0.7 1.2]));

insetB = axes('Units', 'inches', 'Position', [0.45 5.66 5.95 0.94]);

panelC = struct('top_left', axes('Units', 'inches', 'Position', [3.125 4.46 0.7 1.2]), ...
    'top_right', axes('Units', 'inches', 'Position', [4.21 4.46 0.7 1.2]), ...
    'middle_left', axes('Units', 'inches', 'Position', [3.125 2.45 0.7 1.2]), ...
    'middle_right', axes('Units', 'inches', 'Position', [4.21 2.45 0.7 1.2]), ...
    'bottom_left', axes('Units', 'inches', 'Position', [3.125 0.44 0.7 1.2]), ...
    'bottom_right', axes('Units', 'inches', 'Position', [4.21 0.44 0.7 1.2]));

insetC = struct('top', axes('Units', 'inches', 'Position', [3.825 4.46 0.385 1.2]), ...
    'middle', axes('Units', 'inches', 'Position', [3.825 2.45 0.385 1.2]), ...
    'bottom', axes('Units', 'inches', 'Position', [3.825 0.44 0.385 1.2]));

panelD = struct('top', axes('Units', 'inches', 'Position', [5.7 4.46 0.7 1.2]), ...
    'middle', axes('Units', 'inches', 'Position', [5.7 2.45 0.7 1.2]), ...
    'bottom', axes('Units', 'inches', 'Position', [5.7 0.44 0.7 1.2]));

panelE = axes('Units', 'inches', 'Position', [5.25 7 1 1]);

panelF = axes('Units', 'inches', 'Position', [0.52 0.44 1.75 1.2]);


%% Panel A

params = struct('Chr', 3, 'X', [75 115], 'Y', [-1.8 2], 'IRs_shown', [232 240 247 310]);

plot(aggregate_S_G1{params.Chr}(:, 1) ./ 1e6, aggregate_S_G1{params.Chr}(:, 2), '.', ...
    'Color', s_dark, 'Parent', panelA.top)
set(panelA.top, 'XLim', params.X, 'XTick', [], 'YLim', params.Y)
plot(single_cell_IRs{params.Chr}(params.IRs_shown, 4) ./ 1e6 * ones(1, 2), params.Y, 'k', ...
    'LineWidth', 0.7, 'Parent', panelA.top)
ylabel(panelA.top, 'RT')

num_cells = sum(is_included_chr(params.Chr, :));
[Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(params.Chr, :)));
r = replication_state_filtered{params.Chr}(:, is_included_chr(params.Chr, :));
imagesc(genome_windows{params.Chr}(:, 3)./ 1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
    'Parent', panelA.bottom)
set(panelA.bottom, 'YDir', 'reverse', 'XLim', params.X, 'YLim', [1 num_cells], ...
    'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
colormap(panelA.bottom, sc_colormap)
plot(single_cell_IRs{params.Chr}(params.IRs_shown, 4) ./ 1e6 * ones(1, 2), [1 num_cells], ...
    'k', 'LineWidth', 0.7, 'Parent', panelA.bottom)
xlabel(panelA.bottom, ['Chromosome ' num2str(params.Chr) ' Coordinate, Mb'], ...
    'BackgroundColor', 'white')
ylabel(panelA.bottom, '% Replicated')
yyaxis(panelA.bottom, 'right')
set(panelA.bottom, 'YColor', 'k', 'YTick', [])
ylabel(panelA.bottom, [num2str(num_cells) ' cells'])

posA = panelA.bottom.Position(3);
posB = insetB.Position(3);
new_X = (params.X(2) - params.X(1)) / posA * posB;
set(insetB, 'XLim', params.X(1) + [0 new_X], 'YLim', [0 1], 'Visible', 'off')

%% Panels B, C, and D (example IRs)

params = struct('panel', {panelB.top_left panelB.top_right panelB.bottom_left ...
    panelB.bottom_right panelC.top_left panelC.middle_left panelC.bottom_left panelD.top ...
    panelD.middle panelD.bottom}, 'Chr', {3 3 5 9 3 12 5 3 4 10}, ...
    'IR_shown', {232 240 27 125 247 246 15 310 57 131}, 'inset', [], 'inset_Y', NaN, ...
    'AB_X', NaN, 'X', {[78.07 80.07] [80.82 82.82] [8.92 9.92] [75.27 76.27] [83.16 85.16] ...
      [87.3850 87.7850] [3.585 4.985] [109.79 110.79] [19.68 21.68] [52.96 54.16]});

params(1).AB_X = 75;
params(2).AB_X = 86;
params(5).AB_X = 100;
params(8).AB_X = 126;

params(5).inset = insetC.top;
params(6).inset = insetC.middle;
params(7).inset = insetC.bottom;

for p = 1:size(params, 2)
    
    parent = params(p).panel;
    
    X = params(p).X;
    
    index = is_included_chr(params(p).Chr, :);
    w = floor(single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4)/20e3) + 1;
    index = index & ~isnan(replication_state_filtered{params(p).Chr}(w, :));
    num_cells = sum(index);
    [Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(index));
    
    r = replication_state_filtered{params(p).Chr}(:, index);
    imagesc(genome_windows{params(p).Chr}(:, 3)./ 1e6, 1:num_cells, r', 'AlphaData', r', ...
        'Parent', parent)
    set(parent, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 num_cells], 'YTick', Yticks(2:2:end), ...
        'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(parent, sc_colormap)
    plot(single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 .* ones(1,2), ...
        [1 num_cells], 'k', 'LineWidth', 0.7, 'Parent', parent)
    xlabel(parent, ['Chr' num2str(params(p).Chr) ', Mb'])
    
    if ~ismember(p, [2 4])
        ylabel(parent, '% Replicated')
    end
    
    if late_IR_class{params(p).Chr}(params(p).IR_shown, 1)
        title(parent, 'Const. Late');
    elseif late_IR_class{params(p).Chr}(params(p).IR_shown, 2)
        title(parent, 'Early+Rare');
    elseif late_IR_class{params(p).Chr}(params(p).IR_shown, 3)
        title(parent, 'Throughout S');
    end
    
    if ismember(p, [1 2 5 8])
        t = get(parent, 'Title');
        t_pos = t.Position;
        set(t, 'Position', [t_pos(1) -165], 'BackgroundColor', 'white')
    end
    
    if ~isempty(params(p).inset)
        set(params(p).inset, 'YDir', 'reverse', 'XLim', [0 1], 'YLim', [1 num_cells], ...
            'Visible', 'off')
        params(p).inset_Y = num_cells;
    end
    
end

% Zoom in
params(5).panel = panelC.top_right;
params(6).panel = panelC.middle_right;
params(7).panel = panelC.bottom_right;

for p = 5:7
    parent = params(p).panel;

    if p == 6
        X = [87.185 87.985];
    else
        X = single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4)/1e6 + [-0.5 0.5];
    end
    
    index = replication_tracks{params(p).Chr}(:, 8) == params(p).IR_shown;
    barcodes = replication_tracks{params(p).Chr}(index, 1);
    
    if length(barcodes) > 5
        barcodes = barcodes(1:5);
    end
    
    YLabels = percent_replicated_filtered(barcodes) .* 100;

    r_tmp = replication_state_filtered{params(p).Chr}(:, is_included_chr(params(p).Chr, :));

    barcodes_shown = find(is_included_chr(params(p).Chr, :));

    num_cells = length(barcodes);

    pos = get(parent, 'Position');
    pos(4) = 0.24 * num_cells;
    pos(2) = pos(2) + (0.24 * (5 - num_cells));
    set(parent, 'Position', pos)

    col = (num_cells * 11) + (2 * (num_cells-1));
    r = NaN(size(r_tmp, 1), col);
    Y = NaN(num_cells, 1);
    counter = 1;
    for b = 1:num_cells
        barcode = find(barcodes_shown == barcodes(b));
        Y(b) = barcode;
        barcode = barcode +  [-5 5];
        r(:, counter:counter+10) = r_tmp(:, barcode(1):barcode(2));
        
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
    colormap(parent, sc_colormap)
    
    plot(single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 * ones(1,2), [1 63], ...
        'k', 'LineWidth', 0.7, 'Parent', parent)
    
    Y(:, 2) = [0 params(p).inset_Y / 1.2 * pos(4)];
    for l = 1:2
        plot([0 1], Y(l, :), 'k:', 'LineWidth', 1, 'Parent', params(p).inset)
    end

end

%% Connect A and B

for p = [1 2 5 8]
    plot([single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) /1e6 params(p).AB_X], [1 0], ...
        'k', 'LineWidth', 0.7, 'Parent', insetB)
    plot([single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) /1e6 params(p).AB_X + 6.8], ...
        [1 0], 'k', 'LineWidth', 0.7, 'Parent', insetB)
end

uistack(insetB, 'bottom')

%% Panel E

params = struct('text_pos', {[-0.1857 -0.0751], [0.3110 -0.4009], [0.4172 0.066], ...
    [0.2227 1.068]}, ...
    'color', {'#66C2A5', '#E78AC3', '#8DA0CB', '#FC8D62'});

p = pie([sum(flat_p_fired > 0.5) sum(flat_classes, 1)], ones(1,4), 'Parent', panelE);
set(panelE, 'XColor', 'none', 'YColor', 'none')

t = findobj(p, 'Type', 'text');
c = findobj(p, 'Type', 'patch');

for l = 1:4
    t(l).Position = params(l).text_pos;
    t(l).FontName = 'Arial';
    c(l).FaceColor = params(l).color;
end
clearvars t c l

title(panelE, 'IR class', 'Position', [0 1.4631])
legendE = legend(panelE, {'Early in aggregate', 'Throughout S', 'Early + Rare', ...
    'Constitutively Late'});
set(legendE, 'FontSize', 9, 'Units', 'inches', 'Position', [5 6.1875 1.39 0.66])
legendE.ItemTokenSize(1) = 15;

%% Panel F

[~, edges] = histcounts(flat_IRs(flat_p_fired <= 0.5, 11));
for l = [2 3 1]
    h = histogram(flat_IRs(flat_p_fired <= 0.5 & flat_classes(:, l), 11), edges, 'Parent', panelF);
    h.FaceColor = params(l+1).color;
end
set(panelF, 'XDir', 'reverse', 'XLim', [-1.6 0.7], 'YLim', [0 225])

xlabel(panelF, 'Aggregate replication timing')
ylabel(panelF, '# of IRs')
title(panelF, 'Rep. Timing of IR Classes')

%% Annotate panels

params = struct('panel', {panelA.top panelB.top_left panelC.top_left panelD.top panelE panelF}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f'}, ...
    'x', {-0.0746 -0.5545 -0.5545 -0.5545 -0.1806 -0.2302}, ...
    'y', {1.2222 1.0809 1.0809 1.0809 1.25 1.1503});

for p = 1:6
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/Figure5.pdf')
close
