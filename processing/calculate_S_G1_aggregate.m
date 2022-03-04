function aggregate_S_G1 = calculate_S_G1_aggregate(input_files)

    warning('off', 'SPLINES:CHCKXYWP:NaNs')
    load('data/hg37_genome_metadata.mat', 'contigs')

    offset = NaN(length(input_files), 1);
    is_G1 = cell(1, length(input_files));
    for f = 1:length(input_files)
        in = load(['data/intermediate/' input_files{f} '.mat'], 'is_G1');
        is_G1{f} = in.is_G1;
        offset(f) = length(in.is_G1);
    end
    is_G1 = cell2mat(is_G1);
    g1_barcodes = find(is_G1);
    offset = [0; cumsum(offset(1:end-1))];
    reads_per_window = 200;

    g1_windows = cell(22, 1);
    is_masked_window = cell(22, 1);
    raw_counts_20kb = cell(22, 1);

    for Chr = 1:22

        disp(['Chromosome ' num2str(Chr)])
        
        reads = cell(length(input_files), 1);
        for f = 1:length(input_files)
            
            chrom_names = h5info(['data/raw/' input_files{f} '.h5'], '/raw_counts/');
            if any(contains({chrom_names.Datasets.Name},  'chr'))
                input_str = 'chr';
            else
                input_str = '';
            end
            
            reads{f} = double(h5read(['data/raw/' input_files{f} '.h5'], ...
            ['/reads/' input_str num2str(Chr)]));
            reads{f}(:, 1) = reads{f}(:, 1) + offset(f);
        end
        reads = cell2mat(reads);
        reads = sortrows(reads, 2);
        reads(:, 1) = reads(:, 1) + 1;

        % Filter out reads at positions that are not uniquely mappable
        is_map_read = load('data/hg37_mappability_100bp.mat', ['chr' num2str(Chr)]);
        is_map_read = ~is_map_read.(['chr' num2str(Chr)]);
        reads = reads(is_map_read(reads(:, 2)), :);

        % Define G1 windows
        g1_reads = reads(ismember(reads(:, 1), g1_barcodes), 2:3);

        windows = cell(size(contigs{Chr}, 1), 1);
        for contig = 1:size(contigs{Chr},1)
            in_fragment = g1_reads(:,1) > contigs{Chr}(contig,1) & ...
                g1_reads(:,1) <= contigs{Chr}(contig,2);
            if sum(g1_reads(in_fragment, 2)) > reads_per_window
                reads_in_fragment = g1_reads(in_fragment, :);
                reads_in_fragment(:, 3) = 0;
                counter = 0;
                for row = 1:length(reads_in_fragment)
                    counter = counter + reads_in_fragment(row, 2);
                    if counter >= reads_per_window
                        reads_in_fragment(row, 3) = 1;
                        counter = 0;
                    end
                end
                reads_in_fragment(1,3) = 1;
                windows{contig} = find(reads_in_fragment(:, 3));
                windows{contig} = reads_in_fragment(windows{contig}(:, 1), 1);

                windows{contig}(1:end-1, 2) = windows{contig}(2:end, 1) - 1;
                windows{contig}(end, :) = [];

            end
        end

        g1_windows{Chr} = cell2mat(windows);
        g1_windows{Chr}(:, 3) = floor(median(g1_windows{Chr}(:, 1:2), 2));
        g1_windows{Chr} = sortrows(g1_windows{Chr}, 3);

        % Identify noisy windows
        win = g1_windows{Chr}(:, 1:2);
        win = sortrows(win(:));

        half = histcounts(g1_reads(g1_reads(:, 2) == 0.5, 1), win)';
        full = histcounts(g1_reads(g1_reads(:, 2) == 1, 1), win)';
        actual_g1_reads = full(1:2:end) + 0.5 * half(1:2:end);
        is_masked_window{Chr} = actual_g1_reads > reads_per_window + 1 | ...
            actual_g1_reads < reads_per_window -1;

        % Count reads
        s_reads = reads(~ismember(reads(:, 1), g1_barcodes), 2:3);
        s_reads(:, 1) = discretize(s_reads(:, 1), [g1_windows{Chr}(:, 1); g1_windows{Chr}(end, 2)]);
        s_reads = s_reads(~isnan(s_reads(:, 1)), :);

        raw_counts_20kb{Chr} = zeros(size(g1_windows{Chr}, 1), 1);
        for read = 1:size(s_reads, 1)
            raw_counts_20kb{Chr}(s_reads(read, 1), 1) = ...
                raw_counts_20kb{Chr}(s_reads(read, 1), 1) + s_reads(read, 2);
        end

        % Mask noisy windows
        raw_counts_20kb{Chr}(is_masked_window{Chr}, :) = NaN;
    end

    mean_coverage = nanmean(cell2mat(raw_counts_20kb));

    for Chr = 1:22
        raw_counts_20kb{Chr} = raw_counts_20kb{Chr} ./ mean_coverage;
    end

    flat = cell2mat(raw_counts_20kb);
    mu = nanmean(flat);
    sigma = nanstd(flat);

    unsmoothed = cell(22, 1);
    for Chr = 1:22
        unsmoothed{Chr} = (raw_counts_20kb{Chr} - mu) ./ sigma;
        is_extreme = abs(unsmoothed{Chr}) > 3.5;
        unsmoothed{Chr}(is_extreme) = NaN;
    end


    % Smooth the replication timing profile
    smoothed = cell(22, 1);

    for Chr = 1:22

        smoothed{Chr} = NaN(size(unsmoothed{Chr}));

        % Smooth between gaps
        for contig = 1:size(contigs{Chr}, 1)
            index = g1_windows{Chr}(:,1) > contigs{Chr}(contig,1) & ...
                g1_windows{Chr}(:,2) <= contigs{Chr}(contig,2);
            if sum(index) >= 20
                if sum(~isnan(unsmoothed{Chr}(index))) > 2
                    f = csaps(g1_windows{Chr}(index,3), double(unsmoothed{Chr}(index)), 10^-17);
                    smoothed{Chr}(index) = fnval(g1_windows{Chr}(index,3), f);
                end
            end
        end

        index = isnan(unsmoothed{Chr});
        smoothed{Chr}(index) = NaN;
    end

    flat = cell2mat(smoothed);
    mu = nanmean(flat);
    sigma = nanstd(flat);

    for Chr = 1:22
        smoothed{Chr} = (smoothed{Chr} - mu) ./ sigma;
        is_extreme = abs(smoothed{Chr}) > 3.5;
        smoothed{Chr}(is_extreme) = NaN;
    end

    for Chr = 1:22
        d = find([false; diff(g1_windows{Chr}(:, 3)) == 0]);
        while size(d, 1) > 0
            g1_windows{Chr}(d, 3) = g1_windows{Chr}(d, 3) + 0.1;
            d = find([false; diff(g1_windows{Chr}(:, 3)) == 0]);
        end
    end

    aggregate_S_G1 = cell(22, 1);
    for Chr = 1:22
        aggregate_S_G1{Chr} = [g1_windows{Chr}(:, 3) smoothed{Chr}];
    end
    
end

