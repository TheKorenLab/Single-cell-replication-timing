function high_confidence_peaks = call_high_confidence_peaks (aggregate_S_G1)

    gap_length = 1e5;
    min_segment_length = 3e5;
    min_window_per_segment = 30;

    high_confidence_peaks = cell(22, 1);
    for Chr = 1:22

        data = aggregate_S_G1{Chr};
        data(:, 3) = 1:size(data, 1);
        data = data(~isnan(data(:, 2)), :);

        segment_bounds = unique([1; find([0; diff(data(:, 1))] > gap_length); size(data, 1) + 1]);
        segment_replicons = cell(length(segment_bounds) - 1, 1);

        for segment = 1:length(segment_bounds) - 1
            if data(segment_bounds(segment + 1) - 1, 1) - ...
                    data(segment_bounds(segment), 1) >= min_segment_length && ...
                    segment_bounds(segment + 1) - segment_bounds(segment) >= min_window_per_segment

                windows_in_segment = data(segment_bounds(segment, 1): ...
                    segment_bounds(segment + 1, 1) - 1, :);

                peaks = find(islocalmax(windows_in_segment(:, 2)));
                valleys = find(islocalmin(windows_in_segment(:, 2)));

                if isempty(peaks) || isempty(valleys)
                    continue
                end

                if peaks(1) < valleys(1)
                    valleys = [1; valleys]; %#ok<*AGROW>
                end

                if peaks(end) > valleys(end)
                    valleys = [valleys; size(windows_in_segment, 1)];
                end

                replicon_list = [valleys(1:end-1) peaks valleys(2:end)];

                clearvars is_included
                is_included(:, 1) = windows_in_segment(replicon_list(:, 2), 1) - ...
                    windows_in_segment(replicon_list(:, 1), 1) >= 100e3 | ...
                    windows_in_segment(replicon_list(:, 3), 1) - ...
                    windows_in_segment(replicon_list(:, 2), 1) >= 100e3;
                is_included(:, 2) = windows_in_segment(replicon_list(:, 2), 2) - ...
                    windows_in_segment(replicon_list(:, 1), 2) >= 0.05 | ...
                    windows_in_segment(replicon_list(:, 2), 2) - ...
                    windows_in_segment(replicon_list(:, 3), 2) >= 0.05;
                is_included = logical(sum(is_included, 2));

                replicon_list = replicon_list(is_included, 2);
                segment_replicons{segment} = windows_in_segment(replicon_list, 3);
            end
        end

        segment_replicons = cell2mat(segment_replicons);
        high_confidence_peaks{Chr} = aggregate_S_G1{Chr}(segment_replicons, 1);

    end

end
