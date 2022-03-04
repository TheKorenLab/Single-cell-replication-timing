function ratio = calculate_read_count_ratio(counts, states)

    if all(isnan(states))
        ratio = NaN;
        return
    end
    
    % Identify segments of consistent state
    non_nan_values = states;
    non_nan_values(:, 2) = 1:length(non_nan_values);
    non_nan_values = non_nan_values(~isnan(non_nan_values(:, 1)), :);
    
    if all(non_nan_values == 2) || all(non_nan_values == 4)
        ratio = NaN;
        return
    end
    
    is_breakpoint = abs(diff(non_nan_values(:, 1))) == 2;
    segments = [1; non_nan_values(is_breakpoint, 2) + 1];
    segments(:, 2) = [segments(2:end) - 1; length(states)];
    
    % Calculate mean read depth in each segment
    copy_number = NaN(size(segments, 1), 1);
    read_depth = NaN(size(segments, 1), 1);
    for row = 1:size(segments,1)
        
        segment_cn = states(segments(row, 1):segments(row, 2));
        is_assigned_state = ~isnan(segment_cn);
        segment_cn = unique(segment_cn(is_assigned_state));
        assert(length(segment_cn) == 1)
        copy_number(row) = segment_cn;
        
        segment_counts = counts(segments(row,1):segments(row,2));
        segment_counts = segment_counts(is_assigned_state);
        read_depth(row) = nanmean(segment_counts);
    end
    
    ratio = nanmean(read_depth(copy_number == 4)) / nanmean(read_depth(copy_number == 2));
end
