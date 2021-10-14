load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')
load('data/processed/GM12878.mat', 'barcodes')

file_list = dir('data/intermediate/GM12878*.mat');

raw_counts = cell(22, length(file_list));
original_barcodes = cell(length(file_list), 1);
for f = 1:length(file_list)

    load(['data/intermediate/' file_list(f).name], 'raw_counts_20kb', 'g1_windows')
    
    for Chr = 1:22
        d = find([false; diff(g1_windows{Chr}(:, 3)) == 0]);
        while size(d, 1) > 0
            g1_windows{Chr}(d, 3) = g1_windows{Chr}(d, 3) + 0.1;
            d = find([false; diff(g1_windows{Chr}(:, 3)) == 0]);
        end
        
        raw_counts{Chr, f} = interp1(g1_windows{Chr}(:, 3), raw_counts_20kb{Chr}, ...
            genome_windows{Chr}(:, 3));
    end
    
    original_barcodes{f} = 1:size(raw_counts{Chr, f}, 2);
end

sorted_filtered_counts = cell(22, 1);
for Chr = 1:22
    sorted_filtered_counts{Chr} = NaN(size(genome_windows{Chr}, 1), size(barcodes, 1));
    for b = 1:size(barcodes, 1)
        sample = barcodes(b, 1);
        barcode = barcodes(b, 2);
        sorted_filtered_counts{Chr}(:, b) = raw_counts{Chr, sample}(:, barcode);
    end
    
    sorted_filtered_counts{Chr}(~is_mappable{Chr}, :) = NaN;
end

subS_fractions = partition_single_cells(sorted_filtered_counts, [10 100]);

save('-append', 'data/processed/GM12878.mat', 'subS_fractions')
clearvars -except samples
