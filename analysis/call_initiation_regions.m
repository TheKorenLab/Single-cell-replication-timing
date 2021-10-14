function [replication_tracks, single_cell_IRs] = call_initiation_regions(replication_tracks, ...
    aggregate_S_G1, high_confidence_peaks)

    load('data/hg37_genome_metadata.mat', 'genome_windows')

    single_cell_IRs = cell(22, 1);

    for Chr = 1:22

        disp(['Chromosome ' num2str(Chr)])

        % Sort tracks by width and center window
        tracks = replication_tracks{Chr};
        tracks = sortrows(tracks, [6 5]);
        num_tracks = size(tracks, 1);

        % Flag tracks whose boundaries are defined by missing data
        is_bounded_by_gap = tracks(:, 4);
        tracks = tracks(:, [1:3 5:6]);

        % Flag tracks that overlap multiple high-confidence peaks
        contains_multiple_peaks = false(num_tracks, 1);
        overlapped_peak = NaN(num_tracks, 1);
        track_windows = [genome_windows{Chr}(tracks(:, 2), 1) genome_windows{Chr}(tracks(:, 3), 2)];
        for t = 1:num_tracks

            is_overlapped_peak = high_confidence_peaks{Chr} > track_windows(t, 1) & ...
                high_confidence_peaks{Chr} < track_windows(t, 2);

            if sum(is_overlapped_peak) == 1
                overlapped_peak(t) = find(is_overlapped_peak);
            elseif sum(is_overlapped_peak) > 1
                contains_multiple_peaks(t) = true;
            end
        end

        % Add uncertainty around track edges
        blurred_tracks = tracks(:, 2:3) + [-2.5 2.5];
        blurred_tracks(:, 3) = mean([blurred_tracks(:, 1) tracks(:, 3)], 2);
        blurred_tracks(:, 4) = mean([tracks(:, 2) blurred_tracks(:, 2)], 2);

        % Seed the IR list with the aggregate peaks
        putative_IRs = [floor(high_confidence_peaks{Chr}/20e3) + 1 + [-2.5 2.5]; NaN(num_tracks, 2)];
        num_IRs = size(putative_IRs, 1);
        cluster_id = NaN(num_tracks, 1);

        % Assign segments to clusters, prioritizing short segments
        for t = 1:num_tracks

            if is_bounded_by_gap(t) || contains_multiple_peaks(t) % skip these segments
                continue
            end

            % Identify existing clusters that overlap the segment
            left = max([repmat(blurred_tracks(t, 1), num_IRs, 1), putative_IRs(:, 1)], [], 2);
            left(all(isnan(putative_IRs), 2)) = NaN;

            right = min([repmat(blurred_tracks(t, 2), num_IRs, 1), putative_IRs(:, 2)], [], 2);
            right(all(isnan(putative_IRs), 2)) = NaN;

            if tracks(t, 5) <= 10
                min_overlap = 0;
            else
                min_overlap = 1;
            end

            is_overlapping = right - left > min_overlap;
            num_overlaps = sum(is_overlapping);

            if num_overlaps > 1 % segment overlaps multiple clusters
                continue
            end

            if num_overlaps == 0 % seed a new cluser
                new_IR = find(isnan(putative_IRs(:, 1)), 1);
                putative_IRs(new_IR, :) = blurred_tracks(t, 3:4); % center position
                cluster_id(t) = new_IR;
                continue
            end

            % Now we check that the track center overlaps the center of at least one track already
            % assigned to the same cluster

            overlapped_cluster = find(is_overlapping);
            is_in_same_cluster = cluster_id == overlapped_cluster;

            if sum(is_in_same_cluster) == 0 % first segment being added to a seeded cluster
                cluster_id(t) = overlapped_cluster;
                continue
            end

            % Check for overlap

            if any(tracks(is_in_same_cluster, 1) == tracks(t, 1))
                continue
            end

            clustered_tracks = blurred_tracks(is_in_same_cluster, :);
            left = max([repmat(blurred_tracks(t, 3), size(clustered_tracks, 1), 1), ...
                clustered_tracks(:, 3)], [], 2);
            right = min([repmat(blurred_tracks(t, 4), size(clustered_tracks, 1), 1), ...
                clustered_tracks(:, 4)], [], 2);

            if sum(right - left >= 0) == 0 % non-overlapping centers
                continue
            end

            if tracks(t, 5) <= 10
                clustered_tracks(end+1, :) = blurred_tracks(t, :); %#ok<SAGROW>
                putative_IRs(overlapped_cluster, :) = [min(clustered_tracks(:, 3)) ...
                    max(clustered_tracks(:, 4))];
            end

            cluster_id(t) = overlapped_cluster;
        end

        % Prune IR list
        putative_IRs = putative_IRs(all(~isnan(putative_IRs(:, 1:2)), 2), :);
        num_IRs = size(putative_IRs, 1);

        putative_IRs(:, 3) = 1:num_IRs;
        putative_IRs(:, 4) = 0;
        for row = 1:num_IRs
            putative_IRs(row, 4) = sum(cluster_id == row);
        end
        putative_IRs = sortrows(putative_IRs, 1);

        % Re-assign tracks to IRs
        cluster_id = NaN(num_tracks, 1);
        for t = 1:num_tracks

            if is_bounded_by_gap(t) || contains_multiple_peaks(t)
                continue
            end

            left = max([repmat(blurred_tracks(t, 1), size(putative_IRs, 1), 1) ...
                putative_IRs(:, 1)], [], 2);
            right = min([repmat(blurred_tracks(t, 2), size(putative_IRs, 1), 1) ...
                putative_IRs(:, 2)], [], 2);

            is_overlapped_IR = right - left > 0;

            if sum(is_overlapped_IR) ~= 1
                continue
            end

            overlapped_IR = putative_IRs(is_overlapped_IR, 3);

            if overlapped_IR ~= overlapped_peak(t) && ~isnan(overlapped_peak(t))
                continue
            end

            cluster_id(t) = find(is_overlapped_IR);
        end

        for row = 1:num_IRs
            putative_IRs(row, 4) = sum(cluster_id == row);
        end

        % Re-number the clusters
        putative_IRs(:, 3) = 1:size(putative_IRs, 1);
        putative_IRs = putative_IRs(putative_IRs(:, 4) > 0, :);
        num_IRs = size(putative_IRs, 1);

        cluster_id(:, 2) = NaN;
        for row = 1:size(putative_IRs, 1)
            index = cluster_id(:, 1) == putative_IRs(row, 3);
            cluster_id(index, 2) = row;
        end
        cluster_id = cluster_id(:, 2);

        % Convert replication track list to chromosome coordinates
        tracks(:, 2:3) = [genome_windows{Chr}(tracks(:, 2), 1) genome_windows{Chr}(tracks(:, 3), 2)];
        tracks(:, 4) = round(mean(tracks(:, 2:3), 2));
        tracks(:, 5) = tracks(:, 5) .* 20e3;
        tracks(:, 6) = cluster_id;
        tracks(:, 7) = overlapped_peak;

        % Convert IR list to chromosome coordinates
        new_IR_list = NaN(num_IRs, 10);
        new_IR_list(:, 1) = 1:num_IRs;

        for id = 1:num_IRs

            is_in_cluster = cluster_id == id;
            tracks_in_cluster = tracks(is_in_cluster, [4 7]);
            new_IR_list(id, 2:6) = prctile(tracks_in_cluster(:, 1), [0 100 50 25 75]);

            overlapped_peak = unique(tracks_in_cluster(:, 2));
            overlapped_peak = overlapped_peak(~isnan(overlapped_peak));

            if ~isempty(overlapped_peak)
                new_IR_list(id, 9) = overlapped_peak;
            end

            new_IR_list(id, 10) = size(tracks_in_cluster, 1);
        end

        new_IR_list(:, 7) = new_IR_list(:, 6) - new_IR_list(:, 5) + 20e3;

        index = ~isnan(aggregate_S_G1{Chr}(:, 2));
        new_IR_list(:, 11) = interp1(aggregate_S_G1{Chr}(index, 1), aggregate_S_G1{Chr}(index, 2), ...
            new_IR_list(:, 4));

        for id = 1:num_IRs
            dist = abs(new_IR_list(id, 4) - new_IR_list(:, 4));
            dist(id) = NaN;
            new_IR_list(id, 8) = min(dist);
        end

        single_cell_IRs{Chr} = new_IR_list;

        % Save replication tracks
        tracks = tracks(~is_bounded_by_gap, :);
        replication_tracks{Chr} = sortrows(tracks, [1 4]);

    end
    
end
