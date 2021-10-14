load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed/reference_bulk_profiles.mat', 'ref')

data = struct;
for sample = 1:9
    data.(samples{sample}) = load(['data/processed/' samples{sample} '.mat'], ...
        'replication_state_filtered', 'percent_replicated_filtered', 'is_included_chr');
end

cell_line_colors = {'#41AE76','#238B45', '#005824', '#74A9CF', '#2B8CBE', '#045A8D', '#F768A1', ...
    '#C51B8A', '#7A0177'};

%% Figure skeleton

figureS5 = figure;
set(figureS5, 'Position', [25 9 6.5 6.5])

x = linspace(0.45, 5.25, 4);
y = [5.72 4.84 3.96 3.08 2.2 1.32 0.44];
height = [0.42 repmat(0.83, [1 6])];
    
panels = struct;
panel_names = {'bulk', 'a', 'b', 'c', 'd', 'e', 'f'};

for xi = 1:4
    for yi = 1:7
        panels(xi).(panel_names{yi}) = axes('Units', 'inches', ...
            'Position', [x(xi) y(yi) 0.9 height(yi)]);
    end
end

delete(panels(3).f)
delete(panels(4).f)

%% Plot

titles = {'LCL vs. ESC', 'MCF-7 + ESC', 'RKO-Specific'};
sc_colors = {s_light{1} s_light{1} s_light{1} s_light{2} s_light{2} s_light{2} ...
    s_light{3} s_light{3} s_light{3}};
examples = struct('Chr', {1 2 1 4}, 'X', {[194.5 200] [180 184] [79 85] [175 180]}, ...
    'title', {1 1 2 3}, 'samples', {1:6 1:6 [1 4 7:9] [1 4 7:9]}, 'bulk', {[1 2] [1 2] 1:5 1:5}, ...
    'peak', {[196.1362 198.5757] 182.1423 [81.5372 81.9271] 177.8});
sample_names = samples(1:9);
sample_names{7} = 'HCT-116';
sample_names{9} = 'MCF-7';

for e = 1:4

    col = panels(e);

    % bulk
    parent = col.bulk;
    Chr = examples(e).Chr;
    X = examples(e).X;

    for b = 1:length(examples(e).samples)
        bulk = ref.(samples{examples(e).samples(b)});
        c = cell_line_colors{examples(e).samples(b)};
        plot(bulk{Chr}(:, 1) ./ 1e6, bulk{Chr}(:, 2), '.', 'Color', c, 'HandleVisibility', 'off', ...
            'Parent', parent)
    end
    for o = 1:length(examples(e).peak)
        plot(examples(e).peak(o) * ones(1, 2), [-2 2], 'k', 'LineWidth', 0.7, ...
            'HandleVisibility', 'off', 'Parent', parent)
    end

    set(parent, 'XLim', X, 'XTick', [], 'YLim', [-2 2], 'YTick', [-2 0 2])
    ylabel(parent, 'RT')
    title(parent, titles{examples(e).title})

    % cell lines

    for l = 1:length(examples(e).samples)

        parent = col.(panel_names{l+1});
        s_index = examples(e).samples(l);
        c = sc_colors{s_index};

        s = data.(samples{s_index});

        index = s.is_included_chr(Chr, :);
        num_cells = sum(index);
        [Yticks, YLabels] = get_heatmap_yticks(s.percent_replicated_filtered(index));
        r = s.replication_state_filtered{Chr}(:, index);
        imagesc(genome_windows{Chr}(:, 3) ./1e6, 1:num_cells, r', 'AlphaData', ~isnan(r'), ...
            'Parent', parent)
        set(parent, 'YDir', 'reverse', 'XLim', X, 'YLim', [1 num_cells], ...
            'CLim', [2 4], 'YTick', Yticks(2:2:end), 'YTickLabel', YLabels(2:2:end), 'Box', 'off')
        colormap(parent, [convert_hex(g1_light); convert_hex(c)])
        ylabel(parent, '% Rep.')

        if l < length(examples(e).samples)
            set(parent, 'XTick', [])
        else
            xlabel(parent, ['Chr' num2str(Chr) ', Mb'])
        end

        for o = 1:length(examples(e).peak)
            plot(examples(e).peak(o) * ones(1, 2), [1 num_cells], 'k', 'LineWidth', 0.7, ...
                'HandleVisibility', 'off', 'Parent', parent)
        end

        if strcmp(samples{examples(e).samples(l)}, 'HCT116')
            set(parent, 'YTick', Yticks([2 5 8]), 'YTickLabel', YLabels([2 5 8]))
        end
        
        yyaxis(parent, 'right')
        set(parent, 'YColor', 'k', 'YTick', [])
        ylabel(parent, sample_names{s_index})

    end
end

%% Legend

parent = panels(4).bulk;
for b = [1 4 7 2 5 8 3 6 9]
    c = cell_line_colors{b};
    plot(0, 0, 'LineWidth', 2, 'Color', c, 'DisplayName', samples{b}, 'Parent', parent)
end

legendD = legend(panels(4).bulk);
set(legendD, 'Units', 'inches', 'Position', [3.65 0.2 2.5 0.5], 'FontSize', 9, ...
    'Orientation', 'horizontal', 'NumColumns', 3)
legendD.ItemTokenSize(1) = 15;

%% Annotate panels

params = struct('panel', {panels(1).bulk panels(3).bulk}, 'text', {'a', 'b'}, 'x', -0.3721, ...
    'y', 1.4333);

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS5.pdf')
close
