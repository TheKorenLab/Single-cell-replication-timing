load_single_cell_project
load('data/hg37_genome_metadata.mat', 'is_mappable')

file_list = dir('data/intermediate/*.mat');
for f = 1:size(file_list, 1)
    sample_name = strsplit(file_list(f).name, '-');
    file_list(f).sample = sample_name{1};
end

for sample = 1:length(samples)

    disp(samples{sample})
    
    % Combine replicates
    input_files = {file_list(strcmp({file_list.sample}, samples{sample})).name};
    
    replication_state = cell(22, length(input_files));
    barcodes = cell(length(input_files), 1);
    for f = 1:length(input_files)
        r = load(['data/intermediate/' input_files{f}], 'replication_state');
        replication_state(:, f) = r.replication_state;
        
        barcodes{f}(:, 2) = 1:size(r.replication_state{1}, 2);
        barcodes{f}(:, 1) = f;
    end
    
    for Chr = 1:22
        replication_state{Chr, 1} = cell2mat(replication_state(Chr, :));
    end
    replication_state = replication_state(:, 1);
    barcodes = cell2mat(barcodes);
    num_cells = size(barcodes, 1);
    clearvars Chr f input_files r
    
    % Calculate percent replicated
    flat_state = cell2mat(replication_state);
    percent_replicated = sum(flat_state == 4, 1) ./ sum(~isnan(flat_state), 1);
    
    % Consistency of copy number across chromosomes
    mu_chr = NaN(22, num_cells);
    for Chr = 1:22
        mu_chr(Chr, :) = nanmean(replication_state{Chr}, 1);
    end
    sigma_chr = nanstd(mu_chr, 1);
    has_cnvs = any(mu_chr == 2) & any(mu_chr == 4);
    clearvars Chr mu_chr
    
    % Filter and sort
    is_excluded_cell = percent_replicated < 0.03 | isnan(percent_replicated) ...
        | sigma_chr > 0.4 | has_cnvs;
    
    barcodes = barcodes(~is_excluded_cell, :);
    
    percent_replicated_filtered = percent_replicated(~is_excluded_cell);
    [~, sort_order]=sort(percent_replicated_filtered);
    
    replication_state_filtered = cell(22, 1);
    for Chr = 1:22
        replication_state_filtered{Chr} = replication_state{Chr}(:, ~is_excluded_cell);
        replication_state_filtered{Chr} = replication_state_filtered{Chr}(:, sort_order);
    end
    
    barcodes = barcodes(sort_order, :);
   
    % Filter by chromosome, using correlation
    num_cells_analyzed = size(barcodes, 1);
    is_included_chr = false(22, num_cells_analyzed);
    
    r = cell(22, 1);
    parfor Chr = 1:22
        
        r{Chr} = corr(replication_state_filtered{Chr}, 'rows', 'pairwise');
        r{Chr}(r{Chr} == 1) = NaN;
        
        is_negatively_correlated = nanmean(r{Chr}, 1) < 0;
        r{Chr}(is_negatively_correlated, :) = NaN;
        r{Chr}(:, is_negatively_correlated, :) = NaN;
        
        is_excluded = true(num_cells_analyzed, 2);
        
        for barcode = 1:num_cells_analyzed
            
            if sum(replication_state_filtered{Chr}(:, barcode) == 4) == 0 || ...
                    sum(replication_state_filtered{Chr}(:, barcode) == 2) == 0
                continue
            end
            
            neighbors = barcode-5:barcode+5;
            neighbors = neighbors(neighbors > 0 & neighbors <= num_cells_analyzed);
            
            r1 = r{Chr}(neighbors, neighbors);
            index = neighbors == barcode;
            is_excluded(barcode, 1) = sum(sum(r1(index, :) > r1)) == 0;
            
            r2 = r1(:);
            outliers = reshape(isoutlier(r2), size(r1));
            is_excluded(barcode, 2) = (sum(outliers(index, :))) / (size(r1, 1) - 1) > 0.1;
            
        end
        
        is_included_chr(Chr, :) = is_included_chr(Chr, :)' & ...
            ~any([is_negatively_correlated' is_excluded], 2);
        
    end
    
    for Chr = 1:22
        replication_state_filtered{Chr}(:, ~is_included_chr(Chr, :)) = NaN;
    end
    
    num_excluded_chroms = sum(~is_included_chr, 1);
    is_excluded_cell = num_excluded_chroms >= 21;
    barcodes = barcodes(~is_excluded_cell, :);
    
    for Chr = 1:22
        replication_state_filtered{Chr} = replication_state_filtered{Chr}(:, ~is_excluded_cell);
    end
    is_included_chr = is_included_chr(:, ~is_excluded_cell);
    
    % Re-sort by percent replicated
    flat_state = cell2mat(replication_state_filtered);
    percent_replicated_filtered = sum(flat_state == 4, 1) ./ sum(~isnan(flat_state), 1);
    [percent_replicated_filtered, sort_order]=sort(percent_replicated_filtered);
    
    for Chr = 1:22
        replication_state_filtered{Chr} = replication_state_filtered{Chr}(:, sort_order);
    end
    
    is_included_chr = is_included_chr(:, sort_order);
    
    num_cells_analyzed = sum(~is_excluded_cell);
    barcodes = barcodes(sort_order, :);
    
    out_filename = ['data/processed/' samples{sample} '.mat'];
    save('-v7.3', out_filename, 'barcodes', 'replication_state_filtered', ...
        'percent_replicated_filtered', 'is_included_chr', 'num_cells_analyzed')

    clearvars -except file_list samples sample
end
