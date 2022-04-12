cd('~/Desktop/single_cell_manuscript/')
load_single_cell_project

load('data/processed3/GM12878.mat', 'p_cells_fired', 'single_cell_IRs', ...
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
set(figure5, 'Position', [63.5 22.86 18 18.5])

panelA = struct('top', axes('Position', [0.92 17.1 11.3538 1]), ...
    'bottom', axes('Position', [0.92 14.6 11.3538 2.38]));

panelB = struct('top_left', axes(), 'top_right', axes(), ...
    'bottom_left', axes(), 'bottom_right', axes());

insetB = axes('Position', [0.92 12.75 16.58 1.85], 'Visible', 'off');

panelC = struct('top_left', axes(), 'top_right', axes(), ...
    'middle_left', axes(), 'middle_right', axes(), ...
    'bottom_left', axes(), 'bottom_right', axes());

insetC = struct('top', axes(), 'middle', axes(), 'bottom', axes());

panelD = struct('top', axes('Units', 'inches', 'Position', [5.7 4.46 0.7 1.2]), ...
    'middle', axes('Units', 'inches', 'Position', [5.7 2.45 0.7 1.2]), ...
    'bottom', axes('Units', 'inches', 'Position', [5.7 0.44 0.7 1.2]));

panelE = axes('Position', [14.6 15.5 2.33 2.33]);

panelF = axes('Position', [0.92 0.82 5.33 3.1]);

x = [0.92 4.25 8.21 11 15.5];
y = linspace(9.65, 0.82, 3);

pos = cell(5, 3);
for xi = 1:5
    for yi = 1:3
        pos{xi, yi} = [x(xi) y(yi) 2 3.1];
    end
end
pos = cell2mat(pos(:));

panels = {panelB.top_left panelB.top_right panelC.top_left panelC.top_right panelD.top ...
    panelB.bottom_left panelB.bottom_right panelC.middle_left panelC.middle_right panelD.middle ...
    []  [] panelC.bottom_left panelC.bottom_right panelD.bottom};
  
for p = 1:15
    
    if isempty(panels{p})
        continue
    end
    
    set(panels{p}, 'Units', 'centimeters', 'Position', pos(p, :))
end

panels = {insetC.top insetC.middle insetC.bottom};
for p = 1:3
    set(panels{p}, 'Units', 'centimeters', 'Position', [x(3)+2 y(p) 0.79 3.1], 'Visible', 'off');
end

%% Panel A

params = struct('Chr', 3, 'X', [75 115], 'Y', [-1.8 2], 'IRs_shown', [234 242 250 314]);

plot(aggregate_S_G1{params.Chr}(:, 1) ./ 1e6, aggregate_S_G1{params.Chr}(:, 2), '.', ...
    'Color', s_dark, 'Parent', panelA.top)
set(panelA.top, 'XLim', params.X, 'XTick', [], 'YLim', params.Y)
plot(single_cell_IRs{params.Chr}(params.IRs_shown, 4) ./ 1e6 * ones(1, 2), params.Y, 'k', ...
    'LineWidth', 0.6, 'Parent', panelA.top)
ylabel(panelA.top, 'RT')

num_cells = sum(is_included_chr(params.Chr, :));
[Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(params.Chr, :)));
r = replication_state_filtered{params.Chr}(:, is_included_chr(params.Chr, :));
imagesc(genome_windows{params.Chr}(:, 3)./ 1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
    'Parent', panelA.bottom)
set(panelA.bottom, 'YDir', 'reverse', 'XLim', params.X, 'YLim', [0.5 num_cells+0.5], ...
    'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
colormap(panelA.bottom, sc_colormap)
plot(single_cell_IRs{params.Chr}(params.IRs_shown, 4) ./ 1e6 * ones(1, 2), [0.5 num_cells+0.5], ...
    'k', 'LineWidth', 0.6, 'Parent', panelA.bottom)

xlabel(panelA.bottom, ['Chromosome ' num2str(params.Chr) ' Coordinate, Mb'], ...
    'BackgroundColor', 'white')
xl = get(panelA.bottom, 'XLabel');
set(xl, 'Units', 'centimeters');
pos = xl.Position;
pos(2) = -0.5;
set(xl, 'Position', pos)

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
    'IR_shown', {234 242 28 127 250 250 16 314 60 133}, 'inset', [], 'inset_Y', NaN, ...
    'X', {[78.07 80.07] [80.82 82.82] [8.92 9.92] [75.27 76.27] [83.16 85.16] ...
      [87.3850 87.7850] [3.585 4.985] [109.79 110.79] [19.68 21.68] [52.96 54.16]});
  
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
    set(parent, 'YDir', 'reverse', 'XLim', X, 'YLim', [0.5 num_cells+0.5], 'YTick', Yticks(2:2:end), ...
        'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(parent, sc_colormap)
    plot(single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 .* ones(1,2), ...
        [0.5 num_cells+0.5], 'k', 'LineWidth', 0.6, 'Parent', parent)
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
    
    set(parent.Title, 'BackgroundColor', 'white')
    
    if ~isempty(params(p).inset)
        set(params(p).inset, 'YDir', 'reverse', 'XLim', [0 1], 'YLim', [0.5 num_cells+0.5], ...
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
        r(:, counter:counter+10) = r_tmp(:, barcode(1):barcode(2));
        
        if b < size(barcodes, 1)
            r(:, counter+11:counter+12) = NaN;
            counter = counter + 13;
        end
    end
    Y = Y([1 end]);

    imagesc(genome_windows{params(p).Chr}(:, 3) ./ 1e6, 1:col, r', 'AlphaData', ~isnan(r'), ...
        'Parent', parent)
    set(parent, 'XLim', X, 'YLim', [0.5 col+0.5], 'YDir', 'reverse', 'YTick', [6 19 32 45 58], ...
        'YTickLabel', round(YLabels, 0), 'Box', 'off')
    colormap(parent, sc_colormap)
    
    plot(single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4) ./ 1e6 * ones(1,2), [0 63], ...
        'k', 'LineWidth', 0.6, 'Parent', parent)
    
    Y(:, 2) = [0 params(p).inset_Y / 3.1 * pos(4)] + 75;
    for l = 1:2
       plot([0 1], Y(l, :), 'k:', 'LineWidth', 0.6, 'Parent', params(p).inset);
    end

end

%% Connect A and B

X = xlim(insetB);
Mb_per_cm = (X(2)-X(1))/ insetB.Position(3);

panels = {panelB.top_left panelB.top_right [] [] panelC.top_left [] [] panelD.top};

for p = [1 2 5 8]
    
    pos = panels{p}.Position;

    AB_X(1) = X(1) + Mb_per_cm * (pos(1)-insetB.Position(1));
    AB_X(2) = X(1) + Mb_per_cm * (pos(1)+pos(3)-insetB.Position(1));

    for l = 1:2
        plot([single_cell_IRs{params(p).Chr}(params(p).IR_shown, 4)/1e6 AB_X(l)], [1 0], ...
            'k:', 'LineWidth', 0.6, 'Parent', insetB)
    end
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
    t(l).FontSize = 7;
    c(l).FaceColor = params(l).color;
end
clearvars t c l

title(panelE, 'IR class', 'Position', [0 1.35]); 

legendE = legend(panelE, {'Early in aggregate', 'Throughout S', 'Early + Rare', ...
    'Constitutively Late'});
legendE.ItemTokenSize(1) = 15;
set(legendE, 'Units', 'centimeters', 'Position', [14.54 13.64 3 1.32])

%% Panel F

[~, edges] = histcounts(flat_IRs(flat_p_fired <= 0.5, 11));
for l = [2 3 1]
    h = histogram(flat_IRs(flat_p_fired <= 0.5 & flat_classes(:, l), 11), edges, 'Parent', panelF);
    h.FaceColor = params(l+1).color;
end
set(panelF, 'XDir', 'reverse', 'XLim', [-1.6 0.7], 'YLim', [0 225])

xlabel(panelF, 'Aggregate RT')
ylabel(panelF, '# of IRs')
title(panelF, 'RT of IR Classes')

%% Annotate panels

params = struct('panel', {panelA.top panelB.top_left panelC.top_left panelD.top panelE panelF}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f'}, ...
    'x', {-0.0559 -0.2798 -0.2798 -0.2798 0 -0.13}, ...
    'y', {1.2222 1.005 1.005 1.005 1 1.05});

for p = 1:6
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end
   
printFigure('out/v4/Figure5.pdf')
close

clearvars
