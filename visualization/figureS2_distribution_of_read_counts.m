load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/intermediate/GM12878-S-library1-replicate1.mat', 'filtered_counts_20kb', ...
    'replication_state', 'is_G1', 'g1_windows')
replication_state = cell2mat(replication_state);

load('data/processed/GM12878.mat', 'barcodes', 'percent_replicated_filtered')

counts = cell2mat(filtered_counts_20kb);
counts = aggregate_counts(counts, 15);
counts(counts == 0) = NaN;

percent_missing = sum(isnan(counts), 1) ./ size(counts, 1);

g1_barcodes = find(is_G1 & percent_missing < 0.1);
rng(202)
g1_barcodes = randsample(g1_barcodes, 8);
g1_barcodes(2, :) = 0;

s_barcodes(:, 1) = find(~is_G1 & percent_missing < 0.1);
s_barcodes(:, 2) = NaN;
for b = 1:size(s_barcodes, 1)
    index = barcodes(:, 1) == 3 & barcodes(:, 2) == s_barcodes(b);
    assert(sum(index) <= 1)
    
    if sum(index) == 0
        continue
    end
    
    s_barcodes(b, 2) = percent_replicated_filtered(index) * 100;
end
s_barcodes = s_barcodes(~isnan(s_barcodes(:, 2)), :);
s_barcodes = sortrows(s_barcodes, 2, 'ascend');

index = NaN(32, 1);
index(1) = 1;
for r = 2:32
    p = s_barcodes(index(r-1), 2) + 2.6;
    d = abs(s_barcodes(:, 2) - p);
    index(r) = find(d == min(d));
end

s_barcodes = s_barcodes(index, :);

barcodes_shown = [g1_barcodes'; s_barcodes];

chrom_labels = NaN(22, 1);
for Chr = 1:22
    if Chr > 1
        g1_windows{Chr} = g1_windows{Chr} + g1_windows{Chr-1}(end, 2);
        genome_windows{Chr} = genome_windows{Chr} + genome_windows{Chr-1}(end, 2);
    end
    
    chrom_labels(Chr) = nanmedian(genome_windows{Chr}(:, 3));
end
g1_windows = cell2mat(g1_windows);
g1_windows = g1_windows(1:15:end, :);

genome_windows = cell2mat(genome_windows);

clearvars b barcodes cell_line_names Chr d g1_barcodes index is_G1 p percent_missing ...
    percent_replicated_filtered r s_barcodes samples

% Figure skeleton
figureS2 = figure;
set(figureS2, 'Units', 'inches', 'Position', [25 9 6.5 10])

pos_x = linspace(0.355, 5.875, 8);
pos_y = linspace(9.3, 1.3, 5);

pos = cell(8, 4);
for xi = 1:8
    for yi = 1:5
        pos{xi, yi} = [pos_x(xi) pos_y(yi)];
    end
end

pos = cell2mat(reshape(pos, [40 1]));

% Plot histograms
yi = 1;
for s = 1:40
    axes('Units', 'inches', 'Position', [pos(s, :) 0.54 0.54])
    
    [y, edges] = histcounts(counts(:, barcodes_shown(s, 1)));
    x = [edges(1:end-1); edges(2:end)];
    x = mean(x, 1);
    
    y = y ./ sum(y) .* 100;
    bar(x, y, 1, 'FaceColor', '#BDBDBD', 'EdgeColor', '#636363')
    ylim([0 max(y) + 1.5])
    
    if barcodes_shown(s, 2) == 0
        title('G1/G2')
    else
        title([num2str(barcodes_shown(s, 2), '%0.1f') '% rep'])
    end
    
    xlabel('Counts')
    
    if ismember(s, 1:8:40)
        ylabel('% windows')
        yl = get(gca, 'YLabel');
        set(yl, 'Units', 'inches')
        yl_pos = get(yl, 'Position');
        yl_pos(1) = -0.24;
        set(yl, 'Position', yl_pos)
        
        axes('Units', 'inches', 'Position', [pos_x(1) pos_y(yi)-0.72 6.06 0.33])
        plot(g1_windows(:, 3)./1e6, counts(:, barcodes_shown(s, 1)), 'k.', 'MarkerSize', 2)
        ymin = prctile(counts(:, barcodes_shown(s, 1)), 0.5);
        ymax = prctile(counts(:, barcodes_shown(s, 1)), 99.9);
        set(gca, 'XTick', [], 'XLim', [0 g1_windows(end, 3)]/1e6, ...
            'YLim', [ymin - 10 ymax + 15], 'YTick', [50 100])
        ylabel('Counts')
        yl = get(gca, 'YLabel');
        set(yl, 'Units', 'inches')
        yl_pos = get(yl, 'Position');
        yl_pos(1) = -0.24;
        set(yl, 'Position', yl_pos)
        
        axes('Units', 'inches', 'Position', [pos_x(1) pos_y(yi)-0.98 6.06 0.15])
        plot(genome_windows(:, 3)./1e6, replication_state(:, barcodes_shown(s, 1)), 'k.', ...
            'MarkerSize', 2)
        set(gca, 'XLim', [0 genome_windows(end, 3)]/1e6, 'XTick', chrom_labels/1e6, ...
            'XTickLabel', 1:22, 'YTick', [2 4], 'YLim', [1 5])
        xlabel('Chromosome')
        ylabel('CN')
        yl = get(gca, 'YLabel');
        set(yl, 'Units', 'inches')
        yl_pos = get(yl, 'Position');
        yl_pos(1) = -0.24;
        set(yl, 'Position', yl_pos)
        
        yi = yi + 1;
    end
    
end
 
for yi = 1:5
    annotation('arrow', 'Units', 'inches', 'X', [0.855 1], 'Y', [pos_y(yi)-0.05 pos_y(yi)-0.34])
end

printFigure('out/FigureS2.pdf')
close
clearvars
