function [filtered_counts, states] = infer_single_cell_states(counts, breakpoints)

    assert(size(counts, 2) == 1)
    states = NaN(size(counts));
    
    % Filter raw counts
    filtered_counts = filter_single_cell_counts(counts);
        
    % Remove cells with high autocorrelation after filtering out CNVs
    % These likely have residual CNVs that were not successfully removed.
    non_nan_values = filtered_counts(~isnan(filtered_counts));
    
    if isempty(non_nan_values)
        states(:) = NaN;
        return
    end
    
    seq = aggregate_counts(non_nan_values, 15);
    autocorr_400kb = autocorr(seq);
    
    if autocorr_400kb(21) > 0.15
        states(:) = NaN;
        return
    end
    
    % Intialize HMM
    [initial_emissions, initial_transitions, ...
        initial_states] = initialize_single_cell_hmm(filtered_counts);

    if isempty(initial_states)
        states(:) = 2;
        return
    end

    % Train HMM
    [emissions_matrix, transition_matrix, converged] = ...
        train_single_cell_hmm(filtered_counts, initial_emissions, initial_transitions, initial_states);

    if ~converged
        states(:) = NaN;
        return
    end

    % Decode HMM
    states = decode_single_cell_hmm(filtered_counts, transition_matrix, emissions_matrix, ...
        initial_states, breakpoints);
    
    % Remove cells with unrealistic ratio of read counts between states
    read_count_ratio = calculate_read_count_ratio(filtered_counts, states);
    
    if read_count_ratio < 1.5 || read_count_ratio > 2.5
        states(:) = NaN;
    end

end
