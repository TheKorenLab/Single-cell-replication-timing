function viterbi_path = decode_single_cell_hmm(count_data, transition_matrix, emissions_matrix, ...
    initial_states, breakpoints)

num_states = 2;
num_emissions = max(count_data) + 1;

assert(all(size(transition_matrix) == num_states));
assert(size(emissions_matrix, 1) == num_states);
assert(size(emissions_matrix, 2) == num_emissions);

% Convert to log for numerical stability
log_transitions = log(transition_matrix);
log_emissions = log(emissions_matrix);

% Initialize v
first_value = find(~isnan(count_data), 1, 'first');
v = NaN(num_states, 1);
for state = 1:num_states
    v(state, 1) =  sum(initial_states == state) / length(count_data) * ...
        emissions_matrix(state, count_data(first_value) + 1);
end
prior_v = v;

% Determine probability of each state at each position
trellis = NaN(num_states, length(count_data));

for position = 1:length(count_data)
    
    if isnan(count_data(position))
        continue
    end
    
    if(ismember(position, breakpoints))
        prior_v = NaN(num_states, 1);
        for state = 1:num_states
            prior_v(state, 1) =  sum(initial_states == state) / length(count_data) * ...
                emissions_matrix(state, count_data(position) + 1);
        end
    end
        
    for state = 1:num_states
        
        highest_probability = -inf;
        most_likely_state = 0;
        for row = 1:num_states
            prob = prior_v(row) + log_transitions(row, state);
            if prob > highest_probability
                highest_probability = prob;
                most_likely_state = row;
            end
        end
        
        trellis(state, position) = most_likely_state;
        v(state) = log_emissions(state, count_data(position) + 1) + highest_probability;
    end
    
    prior_v = v;
end

% Trace back through the trellis
[~, final_state] = max(v);
last_value = find(~isnan(count_data), 1, 'last');

viterbi_path = NaN(size(count_data));
viterbi_path(last_value) = final_state;

for position = last_value - 1:-1:1
    if isnan(count_data(position) + 1)
        continue
    end
    
    index = position + 1;
    while isnan(count_data(index)) && index < length(viterbi_path)
        index = index + 1;
    end
    
    viterbi_path(position) = trellis(viterbi_path(index), index);
    
end
 
% Convert to copy-number states
viterbi_path(viterbi_path == 2) = 4;
viterbi_path(viterbi_path == 1) = 2;

end
