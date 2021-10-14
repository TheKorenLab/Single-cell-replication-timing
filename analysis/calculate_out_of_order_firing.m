function firing_order_frequency = calculate_out_of_order_firing(single_cell_IRs, ...
    replication_state_filtered)

    load('data/hg37_genome_metadata.mat', 'genome_windows')

    % Pull out replication state at windows containing IR centers
    state_at_IR = cell(22, 1);
    for Chr = 1:22
        
        IR_centers = NaN(size(single_cell_IRs{Chr}, 1), 1);
        for o = 1:size(single_cell_IRs{Chr}, 1)
            IR_centers(o) = find(genome_windows{Chr}(:, 1) < single_cell_IRs{Chr}(o, 4) & ...
                genome_windows{Chr}(:, 2) >= single_cell_IRs{Chr}(o, 4));
        end
        state_at_IR{Chr} = replication_state_filtered{Chr}(IR_centers, :);

    end

    % Sort IRs from earliest to latest timing
    flat_IR_list = cell2mat(single_cell_IRs);
    [~, sort_order] = sortrows(flat_IR_list, 11, 'descend');

    flat_IR_state = cell2mat(state_at_IR);
    flat_IR_state = flat_IR_state(sort_order, :);

    % Remove cells that are missing too much data
    is_missing_IR = sum(~isnan(flat_IR_state), 1) ./ size(flat_IR_state, 1);
    flat_IR_state(:, is_missing_IR < 0.9) = NaN;

    % Predict which IRs have fired
    num_IRs = nansum(flat_IR_state == 4, 1)';

    expected_IR_state = false(size(flat_IR_state));
    for barcode = 1:size(flat_IR_state, 2)
        expected_IR_state(1:num_IRs(barcode), barcode) = true;
    end

    % Calculate frequency of out-of-order firing
    firing_order_frequency = [sum(flat_IR_state == 4 & expected_IR_state, 2) ...
        sum(flat_IR_state == 2 & ~expected_IR_state, 2) ...
        sum(flat_IR_state == 4 & ~expected_IR_state, 2) ...
        sum(flat_IR_state == 2 & expected_IR_state, 2)];

    firing_order_frequency = firing_order_frequency ./ sum(firing_order_frequency, 2) .* 100;
    
end
