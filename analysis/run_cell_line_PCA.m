load_single_cell_project

% Load data
replication_state = cell(22, 9);
pca_labels = cell(9, 1);
counter = 1;
for sample = 1:9
    load(['data/processed/' samples{sample} '.mat'], 'replication_state_filtered')
    replication_state(:, sample) = replication_state_filtered;
    clearvars replication_state_filtered
    
    new_max = counter + size(replication_state{1, sample}, 2) - 1;
    pca_labels(counter:new_max) = {samples{sample}};
    counter = new_max + 1;
end

replication_state = cell2mat(replication_state);

% Remove noisy windows
p_missing_data_per_window = sum(isnan(replication_state), 2) ./ size(replication_state, 2);
replication_state = replication_state(p_missing_data_per_window < 0.1, :);

% Remove noisy cells
p_missing_data_per_cell = sum(isnan(replication_state), 1) ./ size(replication_state, 1);
replication_state = replication_state(:, p_missing_data_per_cell < 0.1);
pca_labels = pca_labels(p_missing_data_per_cell < 0.1);

% Mask remaining missing data
replication_state(isnan(replication_state)) = eps;

% PCA
[~, pcs] = pca(replication_state', 'NumComponents', 2);

save('data/processed/cell_line_pca.mat', 'pcs', 'pca_labels')
