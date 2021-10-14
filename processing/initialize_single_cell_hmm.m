function [emissions, transitions, initial_states] = initialize_single_cell_hmm(count_data)

    % Index NaNs
    missing_values = isnan(count_data);
    non_nan_values = count_data(~missing_values);

    % Aggregate
    aggregation_width = 15;
    
    seq = aggregate_counts(non_nan_values, aggregation_width);

    % Fit two-component Poisson mixture
    mu = fit_two_component_poisson_mixture(seq);
    
    if isempty(mu)
        emissions = [];
        transitions = [];
        initial_states = [];
        return
    end
    
    mu_ratio = mu(2) / mu(1);

    if mu_ratio < 1.5 || mu_ratio > 2.5
        emissions = [];
        transitions = [];
        initial_states = [];
        return
    end
    
    probability_of_state = NaN(length(seq), 2);
    for m = 1:2
        probability_of_state(:, m) = poisspdf(seq, mu(m));
    end
   
    states = NaN(size(non_nan_values, 1), 1);
    for state = 1:length(mu)
        
        index = find(max(probability_of_state, [], 2) == probability_of_state(:, state));        
        index = (index - 1) * aggregation_width + 1;
        index(:, 2) = min(index(:, 1) + aggregation_width - 1, size(non_nan_values, 1));
        
        for row = 1:size(index, 1)
            states(index(row, 1):index(row, 2)) = state;
        end
        
    end
   
    initial_states = NaN(size(missing_values));
    initial_states(~missing_values) = states;
      
    emissions = NaN(2, max(count_data) + 1);
    for state = 1:2
        emissions(state, :) = poisspdf(0:max(count_data), mean(count_data(initial_states == state)));
    end
    
    transitions = NaN(2, 2);
    for state1 = 1:2
        for state2 = 1:2
            transitions(state1, state2) = sum(initial_states(1:end - 1) == state1 & ...
                initial_states(2:end) == state2);
        end
    end
    
    transitions = transitions ./ sum(transitions, 2);
    
end
