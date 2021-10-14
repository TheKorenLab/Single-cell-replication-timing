function filtered_counts = filter_single_cell_counts(count_data)

    % Index NaNs
    missing_values = isnan(count_data);
    non_nan_values = count_data(~missing_values);

    % Aggregate
    aggregation_width = 15;
    seq = aggregate_counts(non_nan_values, aggregation_width);

    % Identify CNVs
    is_not_cnv = true(size(seq, 1), 1);
    while true
        counts_to_filter = seq(is_not_cnv);
        counts_to_filter = fillmissing(counts_to_filter, 'nearest');

        mu = fit_two_component_poisson_mixture(counts_to_filter);

        if isempty(mu)
            is_not_cnv(:) = false;
            break
        end
        
        probability_of_state = NaN(length(counts_to_filter), 2);
        for m = 1:2
            probability_of_state(:, m) = poisspdf(counts_to_filter, mu(m));
        end

        max_log_probability = -max(log(probability_of_state), [], 2);

        changepoints = findchangepts(counts_to_filter, 'MaxNumChanges', 2, 'Statistic', 'linear');
        
        if isempty(changepoints)
            break
        end
        
        changepoints = [1; changepoints; size(counts_to_filter, 1)]; %#ok<AGROW>

        regions = [changepoints(1:end-1) changepoints(2:end)];
        region_width = regions(:, 2) - regions(:, 1) + 1;
        changepoints = regions(region_width == min(region_width), :)';
        

        if median(max_log_probability(changepoints(1):changepoints(2))) > median(max_log_probability) %   3.5
            index = find(is_not_cnv);
            changepoints = index(changepoints);
            is_not_cnv(changepoints(1):changepoints(2)) = false;
        else
            break
        end

    end

    cnv_windows = find(~is_not_cnv);
    
    % Disaggregate windows
    cnv_windows = (cnv_windows - 1) * aggregation_width + 1;
    cnv_windows(:, 2) = min(cnv_windows(:, 1) + aggregation_width - 1, size(non_nan_values, 1));

    counts = non_nan_values;
    for row = 1:length(cnv_windows)
        counts(cnv_windows(row, 1):cnv_windows(row, 2)) = NaN;
    end

    % Add back in masked windows
    filtered_counts = NaN(size(missing_values));
    filtered_counts(~missing_values) = counts;

end
