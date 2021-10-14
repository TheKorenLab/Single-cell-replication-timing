function [IR_range, cumulative_earliest_timing] = calculate_IR_range(single_cell_IRs, ...
    replication_state_filtered, replication_tracks, percent_replicated_filtered)

    load('data/hg37_genome_metadata.mat', 'genome_windows')

    IR_range = cell(22, 1);

    for Chr = 1:22

        % Pull out windows containing IR centers
        IR_centers = NaN(size(single_cell_IRs{Chr}, 1), 1);
        for o = 1:size(single_cell_IRs{Chr}, 1)
            IR_centers(o) = find(genome_windows{Chr}(:, 1) < single_cell_IRs{Chr}(o, 4) & ...
                genome_windows{Chr}(:, 2) >= single_cell_IRs{Chr}(o, 4));
        end

        % Find unreplicated tracks at overlap each IR
        cn2_tracks = list_cn_segments(replication_state_filtered{Chr}, 2);
        overlapped_IRs = false(size(cn2_tracks, 1), length(IR_centers));
        for track = 1:size(cn2_tracks, 1)

            overlapped_IRs(track, :) = IR_centers > cn2_tracks(track, 2) & ...
                IR_centers < cn2_tracks(track, 3);
            cn2_tracks(track, 4) = sum(overlapped_IRs(track, :));
        end

        % Calculate IR range
        IR_range{Chr} = NaN(size(single_cell_IRs{Chr}, 1), 2);
        is_corroborated = false(size(single_cell_IRs{Chr}, 1), 2);
        for row = 1:length(IR_centers)

            % Find the earliest cell with a replication track attributed to this IR
            if single_cell_IRs{Chr}(row, 10) < 2
                continue
            end

            rep_tracks = replication_tracks{Chr}(replication_tracks{Chr}(:, 8) == row, 1);
            assert(length(rep_tracks) >= single_cell_IRs{Chr}(row, 10))

            rep_tracks = percent_replicated_filtered(rep_tracks);

            IR_range{Chr}(row, 1) = rep_tracks(1);
            is_corroborated(row, 1) = rep_tracks(2) - rep_tracks(1) < 0.1;

            % Find the latest cell that has not replicated the center of this IR
            unrep_tracks = cn2_tracks(overlapped_IRs(:, row), :);

            if isempty(unrep_tracks)
                continue
            end

            min_num_IRs = min(unrep_tracks(:, 4));
            unrep_tracks = unrep_tracks(unrep_tracks(:, 4) == min_num_IRs, 1);
            unrep_tracks = percent_replicated_filtered(unrep_tracks);

            IR_range{Chr}(row, 2) = unrep_tracks(end);
            if length(unrep_tracks) > 1
                is_corroborated(row, 2) = unrep_tracks(end) - unrep_tracks(end-1) < 0.1;
            end
        end

        is_imbalanced = IR_range{Chr}(:, 2) < IR_range{Chr}(:, 1);
        IR_range{Chr}(is_imbalanced, :) = NaN;
        IR_range{Chr}(:, 3) = sum(is_corroborated, 2) == 2;

    end

    % Cumulative earliest firing time
    flat_IR_list = cell2mat(IR_range);
    index = logical(flat_IR_list(:, 3));
    index(isnan(flat_IR_list(:, 2))) = false;

    x = 0:0.01:1;
    cumulative_earliest_timing = NaN(length(x), 2);
    for row = 1:length(x)
        cumulative_earliest_timing(row, 2) = sum(flat_IR_list(index, 1) < x(row)) / sum(index);
    end
    cumulative_earliest_timing(:, 1) = x;

end 
