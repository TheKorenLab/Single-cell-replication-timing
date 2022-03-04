function [g1_windows, is_masked_window, raw_counts_20kb] = define_G1_windows(hdf_file, ...
    g1_barcodes, num_cells)

load('data/hg37_genome_metadata.mat', 'contigs')

g1_windows = cell(22, 1);
is_masked_window = cell(22, 1);
raw_counts_20kb = cell(22, 1);

chrom_names = h5info(hdf_file, '/raw_counts/');
if any(contains({chrom_names.Datasets.Name},  'chr'))
    input_str = 'chr';
else
    input_str = '';
end

for Chr = 1:22
    
    reads = double(h5read(hdf_file, ['/reads/' input_str num2str(Chr)]));
    reads = sortrows(reads, 2);
    reads(:, 1) = reads(:, 1) + 1;
    
    % Filter out reads at positions that are not uniquely mappable
    is_map_read = load('TIGER_hg37_100bp/mappabilityMask.mat', ['chr' num2str(Chr)]);
    is_map_read = ~is_map_read.(['chr' num2str(Chr)]);
    reads = reads(is_map_read(reads(:, 2)), :);
    
    % Define G1 windows
    g1_reads = reads(ismember(reads(:, 1), g1_barcodes), 2:3);
    
    windows = cell(size(contigs{Chr}, 1), 1);
    for contig = 1:size(contigs{Chr},1)
        in_fragment = g1_reads(:,1) > contigs{Chr}(contig,1) & ...
            g1_reads(:,1) <= contigs{Chr}(contig,2);
        if sum(g1_reads(in_fragment, 2)) > 200
            reads_in_fragment = g1_reads(in_fragment, :);
            reads_in_fragment(:, 3) = 0;
            counter = 0;
            for row = 1:length(reads_in_fragment)
                counter = counter + reads_in_fragment(row, 2);
                if counter >= 200
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
    is_masked_window{Chr} = actual_g1_reads > 200 + 1 | ...
        actual_g1_reads < 200 -1;
    
    % Count reads
    reads(:, 2) = discretize(reads(:, 2), [g1_windows{Chr}(:, 1); g1_windows{Chr}(end, 2)]);
    reads = reads(~isnan(reads(:, 2)), :);
    
    raw_counts_20kb{Chr} = zeros(size(g1_windows{Chr}, 1), num_cells);
    for read = 1:size(reads, 1)
        raw_counts_20kb{Chr}(reads(read, 2), reads(read, 1)) = ...
            raw_counts_20kb{Chr}(reads(read, 2), reads(read, 1)) + reads(read, 3);
    end
    
    % Mask noisy windows
    raw_counts_20kb{Chr}(is_masked_window{Chr}, :) = NaN;
    
    % Round to match Poisson expectation
    raw_counts_20kb{Chr} = round(raw_counts_20kb{Chr});
end

end
