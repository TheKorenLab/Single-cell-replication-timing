load_single_cell_project
sc_colormap = [convert_hex(g1_light); convert_hex(s_light{1})];

load('data/hg37_genome_metadata.mat', 'genome_windows')

load('data/processed/GM12878.mat', 'is_included_chr', 'percent_replicated_filtered', ...
    'replication_state_filtered', 'window_IR_frequency', 'replication_tracks', ...
    'single_cell_IRs', 'replication_state_masked', 'aggregate_S_G1')

%% Figure skeleton

figure3 = figure;
set(figure3, 'Position', [25 9 6.5 8])

panelA = struct('top', axes('Units', 'inches', 'Position', [0.45 7.05 3.5 0.6]), ...
    'middle', axes('Units', 'inches', 'Position', [0.45 5.9 3.5 1.1]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.45 5.25 3.5 0.6]));

insetA = axes('Units', 'inches', 'Position', [0.45 4.44 3.5 0.81]);

panelB = struct('left', axes('Units', 'inches', 'Position', [0.45 3.34 1.2 1.1]), ...
    'right', axes('Units', 'inches', 'Position', [2.8 3.34 0.9 1.1]));

panelC = struct('left_top', axes('Units', 'inches', 'Position', [0.45 2.24 1.2 0.7]), ...
    'left_middle', axes('Units', 'inches', 'Position', [0.45 1.19 1.2 1]), ...
    'left_bottom', axes('Units', 'inches', 'Position', [0.45 0.44 1.2 0.7]), ...
    'right', axes('Units', 'inches', 'Position', [2.8 0.44 0.9 2]));

panelD = axes('Units', 'inches', 'Position', [5.07 6.64 1.18 0.95]);

panelE = struct('top', axes('Units', 'inches', 'Position', [4.7 4.79 1.55 1.05]), ...
    'bottom', axes('Units', 'inches', 'Position', [4.985 3.49 0.98 1.05]));

insetE = axes('Units', 'inches', 'Position', [4.7 4.54 1.55 0.25]);

panelF = struct('top', axes('Units', 'inches', 'Position', [4.7 1.74 1.55 1.05]), ...
    'bottom', axes('Units', 'inches', 'Position', [4.825 0.44 1.3 1.05]));

insetF = axes('Units', 'inches', 'Position', [4.7 1.49 1.55 0.25]);

%% Panel A

Chr = 2;
X = [50 84];

set(panelA.top, 'XLim', X, 'XTick', [], 'YLim', [-2.25 2.25])
plot(aggregate_S_G1{Chr}(:,1)./1e6, aggregate_S_G1{Chr}(:,2), '.', 'Color', s_dark, ...
    'Parent', panelA.top)
ylabel(panelA.top, 'RT', 'Position', [47.5 0])
title(panelA.top, 'Chromosomal Distribution of Initiation Events')

num_cells = sum(is_included_chr(Chr, :));
[Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(Chr, :)));
imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, ...
    replication_state_filtered{Chr}(:, is_included_chr(Chr, :))', ...
    'AlphaData', ~isnan(replication_state_filtered{Chr}(:, is_included_chr(Chr, :))'), ...
    'Parent', panelA.middle)
set(panelA.middle, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 num_cells], 'CLim', [2 4], ...
    'XTick', [], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
ylabel(panelA.middle, '% Rep.')
colormap(panelA.middle, sc_colormap)
yyaxis(panelA.middle, 'right')
set(panelA.middle, 'YColor', 'k', 'YTick', [], 'YLim', [1 num_cells])
ylabel(panelA.middle, [num2str(num_cells) '  cells'])

plot(genome_windows{Chr}(:,3)./1e6, window_IR_frequency{Chr}, 'k.', 'Parent', panelA.bottom)
set(panelA.bottom, 'XLim', X, 'YLim', [0 105], 'YTick', [])
xlabel(panelA.bottom, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'], 'BackgroundColor', 'white')
ylabel(panelA.bottom, 'Events', 'Position', [47.5 67.5])

%% Panel B

params = struct('panel', {panelB.left panelB.right}, 'IRs_shown', {177:179 236}, ...
    'title', {'3 IRs' '1 IR'}, 'X', {[57.5 59] [78.6 79.5]});

for p = 1:2
    
    parent = params(p).panel;
    
    barcodes_shown = replication_tracks{Chr}(ismember(replication_tracks{Chr}(:, 6), ...
        params(p).IRs_shown), [1 5]);
    barcodes_shown = sortrows(barcodes_shown, 2);
    barcodes_shown = barcodes_shown(:, 1);
    imagesc(genome_windows{Chr}(:,3)./1e6, 1:length(barcodes_shown), ...
        replication_state_masked{Chr}(:, barcodes_shown)', ...
        'AlphaData', ~isnan(replication_state_masked{Chr}(:, barcodes_shown)'), 'Parent', parent)
    set(parent, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [1 length(barcodes_shown)], ...
        'XTick', [], 'YTick', [], 'Box', 'off')
    colormap(parent, sc_colormap)
    ylabel(parent, [num2str(length(barcodes_shown)) ' cells'])
    title(parent, params(p).title)
    
    for o = 1:length(params(p).IRs_shown)
        plot(single_cell_IRs{Chr}(params(p).IRs_shown(o), 4) * ones(1,2) /1e6, ...
            [1 length(barcodes_shown)], 'k-', 'LineWidth', 0.7, 'Parent', parent)
    end
    
end

%% Panel C

params = struct('panel', {panelC.left_top panelC.left_middle panelC.left_bottom panelC.right}, ...
    'IR_shown', 177, 'X', [57.5 59], 'IRs_shown', {177:179});
params(2).IR_shown = 178;
params(3).IR_shown = 179;
params(4).IR_shown = 236;
params(4).X = [78.6 79.5];
params(4).IRs_shown = 236;

for p = 1:4

    parent = params(p).panel;
    barcodes_shown = replication_tracks{Chr}(:, 8) == params(p).IR_shown & ...
        replication_tracks{Chr}(:, 5) < 1e6;
    barcodes_shown = replication_tracks{Chr}(barcodes_shown, [1 5]);
    barcodes_shown = sortrows(barcodes_shown, 2);
    barcodes_shown = barcodes_shown(:, 1);
    
    if p <= 2
        set(parent, 'XTick', [])
    else
        xlabel(parent, ['Chr' num2str(Chr) ' Coordinate, Mb'])
    end
    
    imagesc(genome_windows{Chr}(:,3)./1e6, 1:length(barcodes_shown), ...
        replication_state_filtered{Chr}(:, barcodes_shown)', ...
        'AlphaData', ~isnan(replication_state_filtered{Chr}(:, barcodes_shown)'), 'Parent', parent)
    
    set(parent, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [1 length(barcodes_shown)], ...
        'YTick', [], 'Box', 'off')
    colormap(parent, sc_colormap)
    ylabel(parent, [num2str(length(barcodes_shown)) ' cells'])
    
    for o = 1:length(params(p).IRs_shown)
        plot(single_cell_IRs{Chr}(params(p).IRs_shown(o), 4) / 1e6 * ones(1,2), ...
            [0 length(barcodes_shown)+1], 'k-', 'LineWidth', 0.7, 'Parent', parent)
    end
    
    yyaxis(parent, 'right')
    set(parent, 'YColor', 'k', 'YTick', [])
    ylabel(parent, [num2str(single_cell_IRs{Chr}(params(p).IR_shown, 7)/1e3) 'kb'])

end

%% Connect left insets

set(insetA, 'XLim', X, 'YLim', [0 1], 'Visible', 'off')
uistack(insetA, 'bottom')

panels = {panelA.top panelA.middle panelA.bottom insetA};

lines = [57.5 59 78.6 79.5];
inset_x = [50 61.75 72.75 81.6];

for l = 1:length(lines)
    for p = 1:3
        y = get(panels{p}, 'YLim');
        plot(lines(l) * ones(1, 2), y, 'k-', 'LineWidth', 0.7, 'Parent', panels{p})
    end
    
    plot([lines(l) inset_x(l)], [1 0], 'k-', 'LineWidth', 0.7, 'Parent', insetA)
end

params = struct('top', {panelB.left panelB.right}, 'bottom', {panelC.left_top panelC.right});

for p = 1:2
    set(params(p).top, 'Units', 'normalized')
    set(params(p).bottom, 'Units', 'normalized')
    
    top_lim = params(p).top.Position(2);
    bottom_lim = params(p).bottom.Position(2) + params(p).bottom.Position(4);
    
    x = params(p).top.Position(1) + (0.5 * params(p).top.Position(3));
    y = [top_lim-0.005 bottom_lim+0.005];
    annotation('arrow', x * ones(1,2), y, 'LineWidth', 0.7, 'Color', 'k');
    
    set(params(p).top, 'Units', 'inches')
    set(params(p).bottom, 'Units', 'inches')
end

%% Panel D

flat_IR_list = cell2mat(single_cell_IRs);
IR_width = flat_IR_list(:, [7 10]);

[counts, edges] = histcounts(IR_width(IR_width(:, 2) >= 5, 1)/1e3);
x = [edges(1:end-1); edges(2:end)];
x = mean(x, 1);

counts = counts ./ sum(counts) .* 100;

bar(x, counts, 1, 'FaceColor', '#BDBDBD', 'EdgeColor', '#636363', 'Parent', panelD);
set(panelD, 'XLim', [20 280], 'YLim', [0 16])
xlabel(panelD, 'IR width, kb')
ylabel(panelD, '% of IRs')
title(panelD, 'IR Width')
text(195, 10, 'Median:', 'Parent', panelD, 'FontSize', 11, 'FontName', 'Arial', ...
    'HorizontalAlignment', 'center');
text(195, 6.8, [num2str(median(IR_width(IR_width(:,2) >= 5, 1)/1e3)) 'kb'], ...
   'Parent', panelD, 'FontSize', 10, 'FontName', 'Arial', 'HorizontalAlignment', 'center');

%% Panel E & F

params = struct('panel', {panelE panelF}, 'Chr', {14 10}, 'IR', {96 181}, 'window', 2, ...
    'title', {'Narrowly Localized IR' 'Broadly Localized IR'}, 'inset', {insetE insetF}, ...
    'inset_x', {[51.9 54.5] [68 72]});

for p = 1:2

    top = params(p).panel.top;
    Chr = params(p).Chr;
    X = single_cell_IRs{Chr}(params(p).IR, 2:3) ./ 1e6 + params(p).window .* [-1 1];
    
    w = floor(single_cell_IRs{Chr}(params(p).IR, 4)/20e3) + 1;    
    is_included_chr(Chr, :) = is_included_chr(Chr, :) ...
        & ~isnan(replication_state_filtered{Chr}(w, :));
    
    num_cells = sum(is_included_chr(Chr, :));
    [Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(Chr, :)));
    imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, ...
        replication_state_filtered{Chr}(:, is_included_chr(Chr, :))', ...
        'AlphaData', ~isnan(replication_state_filtered{Chr}(:, is_included_chr(Chr, :))'), ...
        'Parent', top)
    set(top, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 num_cells], 'CLim', [2 4], ...
        'XTick', [], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    ylabel(top, '% Rep.')
    title(top, params(p).title)
    colormap(top, sc_colormap)
    
    yyaxis(top, 'right')
    set(top, 'YColor', 'k', 'YTick', [], 'YLim', [1 num_cells])
    ylabel(top, [num2str(num_cells) '  cells'])

    set(params(p).inset, 'XLim', X, 'YLim', [0 1], 'Visible', 'off')
    
    bottom = params(p).panel.bottom;
    width = single_cell_IRs{Chr}(params(p).IR, 7);
    width = width .* [-1 1] ./ 1e6;
    X = single_cell_IRs{Chr}(params(p).IR, 2:3) ./ 1e6 + width + [-0.2 0.2];

    yyaxis(top, 'left')
    for t = 1:2
        plot(X(t) * ones(1, 2), [1 num_cells], 'k-', 'LineWidth', 0.7, ...
            'Parent', top)
    end
    
    index = replication_tracks{Chr}(:, 6) == params(p).IR;
    barcodes_shown = replication_tracks{Chr}(index, [1 5]);
    barcodes_shown = sortrows(barcodes_shown, 2);
    barcodes_shown = barcodes_shown(:, 1);
    imagesc(genome_windows{Chr}(:,3)./1e6, 1:length(barcodes_shown), ...
        replication_state_masked{Chr}(:, barcodes_shown)', ...
        'AlphaData', ~isnan(replication_state_masked{Chr}(:, barcodes_shown)'), 'Parent', bottom)
    set(bottom, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 length(barcodes_shown)], 'YTick', [], ...
        'Box', 'off')
    colormap(bottom, sc_colormap)
    xlabel(bottom, ['Chr' num2str(Chr) ', Mb'])
    ylabel(bottom, [num2str(length(barcodes_shown)) ' cells'])
    
    plot(single_cell_IRs{Chr}(params(p).IR, 4) ./ 1e6 * ones(1, 2), [1 num_cells], ...
            'k--', 'LineWidth', 0.7, 'Parent', top)
    for l = 5:6
        plot(single_cell_IRs{Chr}(params(p).IR, l) ./ 1e6 * ones(1, 2), [1 num_cells], ...
            'k--', 'LineWidth', 0.7, 'Parent', bottom) 
    end

    yyaxis(bottom, 'right')
    set(bottom, 'YColor', 'k', 'YTick', [])
    ylabel(bottom, [num2str(single_cell_IRs{Chr}(params(p).IR, 7)/1e3) 'kb'])

    for t = 1:2
        plot([X(t) params(p).inset_x(t)], [1 0], 'k-', 'LineWidth', 0.7, ...
            'Parent', params(p).inset)
    end
    
end

%% Annotate panels
    
params = struct('panel', {panelA.top panelB.left panelC.left_top panelD panelE.top panelF.top}, ...
    'text', {'a', 'b', 'c', 'd', 'e', 'f'}, 'x', {-0.0913 -0.2659 -0.2775 -0.496 -0.1345 -0.1345}, ...
    'y', {1.1839 1.0506 1.06 1.1707 1.1102 1.1102});

for p = 1:6
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/Figure3.pdf')
close
