load_single_cell_project

sc_colormap = [convert_hex(g1_light); convert_hex(s_light{1})];

load('data/hg37_genome_metadata.mat', 'genome_windows')

load('data/processed/GM12878.mat', 'is_included_chr', 'percent_replicated_filtered', ...
    'replication_state_filtered', 'replication_tracks', 'single_cell_IRs', ...
    'replication_state_masked')

figureS2 = figure;
set(figureS2, 'Position', [25 9 6.5 8.94])

x = linspace(0.45, 4.83, 3);
y = [7.59 4.64 1.69];

panel_names = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'};

panels = struct;
p = 0;
for yi = 1:3
    for xi = 1:3
        p = p + 1;
        panels.(panel_names{p}).top = axes('Units', 'inches', 'Position', [x(xi) y(yi) 1.42 1]);
        panels.(panel_names{p}).bottom = axes('Units', 'inches', ...
            'Position', [x(xi)+0.075 y(yi)-1.25 1.3 1]);
        panels.(panel_names{p}).inset = axes('Units', 'inches', ...
            'Position', [x(xi) y(yi)-0.25 1.45 0.25]);
    end
end

params = struct('panel', panel_names, 'Chr', {1 2 3 5 6 7 8 4 1}, ...
    'IR', {331 396 286 21 245 364 400 100 198}, 'window', 2);

for p = 1:size(params, 2)

    top = panels.(params(p).panel).top;
    Chr = params(p).Chr;
    X = single_cell_IRs{Chr}(params(p).IR, 2:3) ./ 1e6 + params(p).window .* [-1 1];
    
    num_cells = sum(is_included_chr(Chr, :));
    [Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(is_included_chr(Chr, :)));
    imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, ...
        replication_state_filtered{Chr}(:, is_included_chr(Chr, :))', ...
        'AlphaData', ~isnan(replication_state_filtered{Chr}(:, is_included_chr(Chr, :))'), ...
        'Parent', top)
    set(top, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 num_cells], 'CLim', [2 4], ...
        'XTick', [], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    ylabel(top, '% Rep.')
    colormap(top, sc_colormap)
    
    yyaxis(top, 'right')
    set(top, 'YColor', 'k', 'YTick', [], 'YLim', [1 num_cells])
    ylabel(top, [num2str(num_cells) '  cells'])

    set(panels.(params(p).panel).inset, 'XLim', X, 'YLim', [0 1], 'Visible', 'off')
    inset_x = X;

    bottom = panels.(params(p).panel).bottom;
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

    offset = get(panels.(params(p).panel).inset, 'XLim');
    width_inset = panels.(params(p).panel).inset.Position(3);
    width_bottom = panels.(params(p).panel).bottom.Position(3);
    offset = (offset(2) - offset(1)) / width_inset * (width_inset - width_bottom) / 2;
    
    inset_x = inset_x + [offset -offset];
    for t = 1:2
        plot([X(t) inset_x(t)], [1 0], 'k-', 'LineWidth', 0.7, ...
            'Parent', panels.(params(p).panel).inset)
    end
    
end

printFigure('out/FigureS2.pdf')
close
