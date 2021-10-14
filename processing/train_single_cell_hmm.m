function [emissions_matrix, transition_matrix, converged] = train_single_cell_hmm(count_data, ...
                                                             initial_emissions, ...
                                                             initial_transitions, initial_states)
    tolerance = 1e-6;
    maxiter = 500;

    num_states = 2;
    num_emissions = max(count_data) + 1;

    assert(all(size(initial_transitions) == num_states));
    assert(size(initial_emissions, 1) == num_states);
    assert(size(initial_emissions, 2) == num_emissions);

    % initialize the counters
    transition_matrix = initial_transitions;
    emissions_matrix = initial_emissions;

    converged = false;
    log_likelihood = 1;

    for iteration = 1:maxiter
        previous_iter = struct('log_likelihood', log_likelihood, ...
            'emissions', emissions_matrix, 'transitions', transition_matrix);
        log_likelihood = 0;

        % Calculate posterior probabilities by Baum-Welch algorithm
        [log_probability_seq, alpha, beta, scalar] = single_cell_baum_welch(count_data, ...
            transition_matrix, emissions_matrix, initial_states);

        log_likelihood = log_likelihood + log_probability_seq;

        % Convert to log for numerical stability
        log_alpha = log(alpha);
        log_beta = log(beta);
        log_transitions = log(transition_matrix);
        log_emissions = log(emissions_matrix);

        % Update transition matrix
        last_position = find(~isnan(count_data), 1, 'last');
        transition_matrix = zeros(size(transition_matrix));
        for state1 = 1:num_states
            for state2 = 1:num_states
                for position = 1:last_position - 1

                    if isnan(count_data(position))
                        continue
                    end

                    index = position + 1;
                    while isnan(count_data(index)) && index < last_position
                        index = index + 1;
                    end

                    t = exp(log_alpha(state1, position) + log_transitions(state1, state2) + ...
                        log_emissions(state2, count_data(index) + 1) + log_beta(state2, index));
                    t = t ./ scalar(index);
                    transition_matrix(state1, state2) = transition_matrix(state1, state2) + t;
                end
            end
        end

        transition_matrix = transition_matrix ./ sum(transition_matrix, 2);

        % Update emissions matrix
        emissions_matrix = zeros(size(emissions_matrix));
        for state = 1:num_states
            for value = 0:max(count_data)
                position = find(count_data == value);
                e = sum(exp(log_alpha(state, position) + log_beta(state, position)));
                emissions_matrix(state, value + 1) = emissions_matrix(state, value + 1) + e;
            end
        end

        emissions_matrix = emissions_matrix ./ sum(emissions_matrix, 2);

        % Handle states with zero transitions
        has_no_transitions = sum(transition_matrix, 2) == 0;
        if any(has_no_transitions)
            transition_matrix(has_no_transitions, :) = 0;
            transition_matrix(sub2ind(size(transition_matrix), ...
                has_no_transitions, has_no_transitions)) = 1;
        end

        % Remove NaNs
        transition_matrix(isnan(transition_matrix)) = 0;
        emissions_matrix(isnan(emissions_matrix)) = 0;

        % Check for convergence
        likelihood_converged = (abs(log_likelihood - previous_iter.log_likelihood) / ...
            (1 + abs(previous_iter.log_likelihood))) < tolerance;
        transitions_converged = norm(transition_matrix - previous_iter.transitions, inf) / ...
            num_states < tolerance;
        emissions_converged = norm(emissions_matrix - previous_iter.emissions, inf) / ...
            num_emissions < tolerance;

       if likelihood_converged && transitions_converged && emissions_converged
            converged = true;
            break
        end
    end
    
end
