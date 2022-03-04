load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed/reference_bulk_profiles.mat', 'ref', 'r_bulk')

% Load single cell data
data = struct;
for sample = 1:9
    data.(samples{sample}) = load(['data/processed/' samples{sample} '.mat'], ...
        'replication_state_filtered', 'percent_replicated_filtered', 'is_included_chr', ...
        'aggregate_S_G1');
end

%% Figure skeleton

figureS7 = figure;
set(figureS7, 'Units', 'inches', 'Position', [25 9 6.5 8.7])

x = [0.45 3.75];

y = [7.94 5.79 3.64 1.49];

panel_names = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};

panels = struct;
p = 0;
for xi = 1:2
    for yi = 1:4
        p = p + 1;
        panels.(panel_names{p}).top = axes('Units', 'inches', 'Position', [x(xi) y(yi) 2.5 0.4]);
        panels.(panel_names{p}).bottom = axes('Units', 'inches', ...
            'Position', [x(xi) y(yi)-1.05 2.5 1]);
    end
end

%% Plot

params = struct('Chr', 2, 'X', [100 150], 'panel', panel_names, 'sample', {2 3 4 5 6 7 8 9}, ...
    'color', {s_light{1} s_light{1} s_light{2} s_light{2} s_light{2} ...
    s_light{3} s_light{3} s_light{3}});

for p = 1:size(params, 2)
    
    top = panels.(params(p).panel).top;
    bottom = panels.(params(p).panel).bottom;
    Chr = params(p).Chr;
    
    bulk = ref.(samples{params(p).sample});
    aggregate_G1S = data.(samples{params(p).sample}).aggregate_S_G1;
    
    is_included_chr = data.(samples{params(p).sample}).is_included_chr;
    percent_replicated_filtered = data.(samples{params(p).sample}).percent_replicated_filtered;
    replication_state_filtered = data.(samples{params(p).sample}).replication_state_filtered;
    
    plot(bulk{Chr}(:, 1) ./ 1e6, bulk{Chr}(:, 2), 'k.', 'Parent', top)
    plot(aggregate_G1S{Chr}(:, 1) ./ 1e6, aggregate_G1S{Chr}(:, 2), '.', ...
        'Color', params(p).color, 'Parent', top)
    set(top, 'XLim', params(p).X, 'XTick', [], 'YTick', [-2 0 2])
    ylabel(top, 'RT')
    title(top, cell_line_names{params(p).sample})
    
    yyaxis(top, 'right')
    set(top, 'YColor', 'k', 'YTick', [])
    ylabel(top, ['r = ' num2str(r_bulk.(samples{params(p).sample}), '%0.2f')])
    

    index = is_included_chr(Chr, :);
    num_cells = sum(index);
    [Yticks, YLabels] = get_heatmap_yticks(percent_replicated_filtered(index));
    r = replication_state_filtered{Chr}(:, index);
    imagesc(genome_windows{Chr}(:, 3) ./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
        'Parent', bottom)
    set(bottom, 'YDir', 'reverse', 'XLim', params(p).X, 'YLim', [0.5 num_cells+0.5], ...
        'CLim', [2 4], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
    colormap(bottom, [convert_hex(g1_light); convert_hex(params(p).color)])
    xlabel(bottom, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    ylabel(bottom, '% Replicated')
    
    if p == 6
        set(bottom, 'YTick', Yticks([2 5 8]), 'YTickLabel', YLabels([2 5 8]))
    end
    
    yyaxis(bottom, 'right')
    set(bottom, 'YColor', 'k', 'YTick', [])
    ylabel(bottom, [num2str(num_cells) ' cells'])
    
end

printFigure('out/FigureS7.pdf')
close
clearvars
