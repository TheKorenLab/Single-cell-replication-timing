load_single_cell_project
load('data/hg37_genome_metadata.mat', 'genome_windows')

file_list = dir('data/intermediate/*.mat');
for f = 1:size(file_list, 1)
    sample_name = strsplit(file_list(f).name, '-');
    file_list(f).sample = sample_name{1};
end

for sample = [1 10]
    input_files = {file_list(strcmp({file_list.sample}, samples{sample})).name};

    load(['data/processed/' samples{sample} '.mat'], 'barcodes')

    filtered_cells = cell(22, length(input_files));
    percent_replicated = cell(1, length(input_files));
    n = NaN(length(input_files), 5);

    for f = 1:length(input_files)

        load(['data/intermediate/' input_files{f}], 'replication_state')

        b = barcodes(barcodes(:, 1) == f, 2);

        n(f, 1) = size(replication_state{1}, 2);
        n(f, 2) = length(b);

        index = 1:size(replication_state{1}, 2);
        index = ~ismember(index, b);

        for Chr = 1:22
            replication_state{Chr} = replication_state{Chr}(:, index);
        end

        flat_states = cell2mat(replication_state);
        p = sum(flat_states == 4, 1) ./ sum(~isnan(flat_states), 1);

        for Chr = 1:22
            filtered_cells{Chr, f} = replication_state{Chr}(:, p > 0);
        end    

        n(f, 3:5) = [sum(p == 0) sum(isnan(p)) sum(p > 0)];
        assert(n(f, 1) == sum(n(f, 2:end)))

        percent_replicated{f} = p(p > 0);

        if sample == 1 && f == 3
            load(['data/intermediate/' input_files{f}], 'raw_counts_20kb', 'g1_windows')
            noisy_cells = cell(22, 1);
            for Chr = 1:22
                noisy_cells{Chr} = raw_counts_20kb{Chr}(:, index);
                noisy_cells{Chr} = noisy_cells{Chr}(:, isnan(p));
            end
            clearvars raw_counts_20kb
        end

        clearvars replication_state
    end

    n = sum(n, 1);

    disp(samples{sample})
    n_labels = {'total', 'analyzed', 'G1/G2', 'noisy', 'removed by QC'};

    disp([num2str(n(1)) ' cells ' n_labels{1}])
    for l = 2:5
        disp([num2str(n(l)) ' (' num2str(100 * n(l)/n(1), '%0.1f') '%) cells ' n_labels{l}])
    end

    if sample == 1
        
        percent_replicated = cell2mat(percent_replicated);
        [~, sort_order] = sort(percent_replicated, 'ascend');
        
        for Chr = 1:22
            p = cell2mat(filtered_cells(Chr, :));
            filtered_cells{Chr, 1} = p(:, sort_order);
        end
        
        filtered_cells = filtered_cells(:, 1);
        
        
        noisy_cells = cell2mat(noisy_cells);
        
        chrom_ticks = NaN(22, 1);
        for Chr = 1:22
            if Chr > 1
                g1_windows{Chr} = g1_windows{Chr} + g1_windows{Chr-1}(end, 2);
            end
            chrom_ticks(Chr, :) = mean(g1_windows{Chr}(:, 3));
        end
        g1_windows = cell2mat(g1_windows);
        
        % Figure
        
        figureR4 = figure;
        set(figureR4, 'Units', 'inches', 'Position', [25 9 6.5 4.3])
        
        y = linspace(3.4861, 1.9861, 3);
        panelA = struct('top', axes('Units', 'inches', 'Position', [0.4 y(1) 5.9 0.6]), ...
            'middle', axes('Units', 'inches', 'Position', [0.4 y(2) 5.9 0.6]), ...
            'bottom', axes('Units', 'inches', 'Position', [0.4 y(3) 5.9 0.6]));
        
        panelB = axes('Units', 'inches', 'Position', [0.4 0.32 5.9 1.1], 'Box', 'off');
        
        % panelA
        params = struct('panel', {panelA.top panelA.middle panelA.bottom}, 'barcode', {5 8 9});
        
        for p = 1:3
            parent = params(p).panel;
            
            counts = aggregate_counts(noisy_cells(:, params(p).barcode), 15);
            ymax = prctile(counts, 99.9);
            
            plot(g1_windows(1:15:end, 3)./1e6, counts, 'k.', 'MarkerSize', 3, 'Parent', parent)
            set(parent, 'XLim', [0 g1_windows(end, 2)]./1e6, 'YLim', [0 ymax])
            
            if p == 3
                set(parent, 'XTickLabels', [1:18 20 22], 'XTick', chrom_ticks([1:18 20 22])/1e6)
                xlabel(parent, 'Chromosome')
            else
                set(parent, 'XTick', [])
            end
            
            ylabel(parent, 'Counts')
        end
        
        % panelB
        Chr = 1;
        X = [0 genome_windows{Chr}(end, 2)/1e6];
        num_cells = size(filtered_cells{Chr}, 2);
        
        imagesc(genome_windows{Chr}(:, 3)./1e6, 1:num_cells, filtered_cells{Chr}', ...
            'AlphaData', ~isnan(filtered_cells{Chr}'), 'Parent', panelB)
        set(panelB, 'XLim', X, 'YLim', [0.5 num_cells+0.5], 'YDir', 'reverse')
        colormap(panelB, [convert_hex(g1_light); convert_hex(s_light{1})])
        xlabel(panelB, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
        ylabel(panelB, 'Cell #')
        
        % Annotate
        params = struct('panel', {panelA.top panelB}, ...
            'text', {'a', 'b'}, 'x', -0.0377, 'y', {1.1395 1.0985});
        
        for p = 1:2
            text(params(p).x, params(p).y, params(p).text, 'Parent', params(p).panel, ...
                'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
        end
        
        printFigure('out/FigureR4.pdf')
        close all

    end
end

clearvars
