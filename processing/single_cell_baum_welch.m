function [log_probability_seq, alpha, beta, scalar] = single_cell_baum_welch(seq, ...
    transition_matrix, emissions_matrix, initial_states)

num_states = size(transition_matrix, 1);
num_emissions = max(seq) + 1;

assert(all(size(transition_matrix) == num_states));
assert(size(emissions_matrix, 1) == num_states)
assert(size(emissions_matrix, 2) == num_emissions);

% Forward procedure
alpha = NaN(num_states, length(seq));
scalar = NaN(1, length(seq));

first_value = find(~isnan(seq), 1, 'first');
for state = 1:num_states
    alpha(state, first_value) = sum(initial_states == state) / length(seq) * ...
        emissions_matrix(state, seq(first_value) + 1);    
end
scalar(first_value) = sum(alpha(:, first_value));
alpha(:, first_value) = alpha(:, first_value) ./ scalar(first_value);

for position = first_value + 1:length(seq)
    if isnan(seq(position))
        continue
    end
    
    index = position - 1;
    while isnan(seq(index)) && index > first_value
        index = index - 1;
    end
    
    for state = 1:num_states
        alpha(state, position) = emissions_matrix(state, seq(position) + 1) .* ...
            (sum(alpha(:, index) .* transition_matrix(:, state)));
    end
    
    scalar(position) = sum(alpha(:, position));
    alpha(:, position) = alpha(:, position) ./ scalar(position);
end

% Backward procedure
last_value = find(~isnan(seq), 1, 'last');
beta = NaN(num_states, length(seq));

beta(:, last_value) = 1;

for position = last_value - 1:-1:1
    if isnan(seq(position))
        continue
    end
    
    index = position + 1;
    while isnan(seq(index)) && index < last_value
        index = index + 1;
    end
    
    for state = 1:num_states
        beta(state, position) = sum(transition_matrix(state, :)' .* beta(:, index) ...
            .* emissions_matrix(:, seq(index) + 1));
        beta(state, position) = beta(state, position) / scalar(index);
    end
        
end

log_probability_seq = nansum(log(scalar));

% 
% seq = [numSymbols+1 seq];
% L = length(seq);
% 
% % This is what we'd like to do but it is numerically unstable
% % warnState = warning('off');
% % logTR = log(tr);
% % logE = log(e);
% % warning(warnState);
% % f = zeros(numStates,L);
% % f(1,1) = 1;
% % % for count = 2:L
% %     for state = 1:numStates
% %         f(state,count) = logE(state,seq(count)) + log(sum( exp(f(:,count-1) + logTR(:,state))));
% %     end
% % end
% % f = exp(f);
% 
% % so we introduce a scaling factor
% fs = zeros(numStates,L);
% fs(1,1) = 1;  % assume that we start in state 1.
% s = zeros(1,L);
% s(1) = 1;
% for count = 2:L
%     if ~isnan(seq(count))
%         c = count - 1;
%         while isnan(seq(c))
%             c = c - 1;
%         end
%         for state = 1:numStates
%             
%             
%             fs(state,count) = e(state,seq(count) + 1) .* (sum(fs(:,c) .*tr(:,state)));
%         end
%     end
%     
%     % scale factor normalizes sum(fs,count) to be 1.
%     s(count) =  sum(fs(:,count));
%     fs(:,count) =  fs(:,count)./s(count);
% end
% s(s == 0) = NaN;
% 
% %  The  actual forward and  probabilities can be recovered by using
% %   f = fs.*repmat(cumprod(s),size(fs,1),1);
% 
% 
% % This is what we'd like to do but it is numerically unstable
% % b = zeros(numStates,L);
% % for count = L-1:-1:1
% %     for state = 1:numStates
% %         b(state,count) = log(sum(exp(logTR(state,:)' + logE(:,seq(count+1)) + b(:,count+1)  )));
% %     end
% % end
% 
% % so once again use the scale factor
% bs = NaN(numStates,L);
% bs(:, end) = 1;
% 
% for count = L-1:-1:1
%     
%     if isnan(seq(count))
%         continue
%     end
%     
%     c = count + 1;
%     if isnan(seq(count + 1))
%         while isnan(seq(c))
%             c = c + 1;
%         end
%     end
%         
%     for state = 1:numStates
%         bs(state,count) = (1/s(c)) * sum( tr(state,:)'.* bs(:,c) .* e(:,seq(c) +1));
%     end
% end
% %  The  actual backward and  probabilities can be recovered by using
% %  scales = cumprod(s, 'reverse'); 
% %  b = bs.*repmat([scales(2:end), 1],size(bs,1),1);
% 
% pSeq = nansum(log(s));
% pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
% pStates(:,1) = [];


% end
