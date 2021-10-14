function assign_replication_states(sample_name)

    load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')
    load(['data/intermediate/' sample_name '.mat'], 'g1_windows', 'is_masked_window', ...
        'raw_counts_20kb')    
    
    % Keep track of chromosome boundaries
    chrom_boundaries = NaN(22, 1);
    for Chr = 1:22
        chrom_boundaries(Chr) = size(g1_windows{Chr}, 1);
    end
    chrom_boundaries = cumsum(chrom_boundaries);
    chrom_boundaries = [1; chrom_boundaries(1:end-1) + 1];

    % Infer replication state by HMM
    flat_raw_counts = cell2mat(raw_counts_20kb);
    num_cells = size(flat_raw_counts, 2);
    flat_masked_windows = cell2mat(is_masked_window);

    flat_filtered_counts = NaN(size(flat_raw_counts));
    flat_replication_inference = NaN(size(flat_raw_counts));

    parfor barcode = 1:num_cells
        disp(['Cell ' num2str(barcode)])
        [flat_filtered_counts(:, barcode), flat_replication_inference(:, barcode)] = ...
            infer_single_cell_states(flat_raw_counts(:, barcode), chrom_boundaries); %#ok<PFOUS>
   
    end

    flat_replication_inference(flat_masked_windows, :) = NaN;

    % Un-flatten chromosomes
    chrom_boundaries(:, 2) = [chrom_boundaries(2:end, 1) - 1; size(flat_replication_inference, 1)];

    filtered_counts_20kb = cell(22, 1);
    replication_inference_g1_windows = cell(22, 1);
    for Chr = 1:22
        index = chrom_boundaries(Chr, :);
        filtered_counts_20kb{Chr} = flat_filtered_counts(index(1):index(2), :);
        replication_inference_g1_windows{Chr} = flat_replication_inference(index(1):index(2), :);
    end

    % Interpolate to fixed-sized windows
    replication_state = cell(22, 1);
    for Chr = 1:22
        replication_state{Chr} = interp1(g1_windows{Chr}(~is_masked_window{Chr}, 3), ...
            replication_inference_g1_windows{Chr}(~is_masked_window{Chr}, :), ...
            genome_windows{Chr}(:, 3));
        replication_state{Chr}(~is_mappable{Chr}, :) = NaN;
        replication_state{Chr}(~ismember(replication_state{Chr}, 2:4)) = NaN;
    end

    % Fill across gaps if the gap is smaller than 100Kb AND the state is the same on both sides
    
    for Chr = 1:22
        
        m = find(~is_mappable{Chr});
        index = [1; find(diff(m) > 1) + 1];
        index(:, 2) = [index(2:end) - 1; size(m, 1)];
        is_missing = [m(index(:, 1)) m(index(:, 2))];
        is_missing(:, 3) = is_missing(:, 2) - is_missing(:, 1) + 1;
        is_missing = is_missing(is_missing(:, 3) <= 5, 1:2);
        
        if isempty(is_missing)
            break
        end
        
        for row = 1:size(is_missing, 1)
            range = is_missing(row, 1)-1:is_missing(row, 2)+1;
            
            if any(range <= 0 | range >= size(is_mappable{Chr}, 1))
                continue
            end
            
            values = replication_state{Chr}(range, :);
            to_fill = values(1, :) == values(end, :);
            for inner = 2:size(values, 1) - 1
                values(inner, to_fill) = values(1, to_fill);
            end
            
            replication_state{Chr}(range, :) = values;
        end
        
    end
    
    save('-append', ['data/intermediate/' sample_name '.mat'], 'filtered_counts_20kb', ...
        'replication_inference_g1_windows', 'replication_state')

end
