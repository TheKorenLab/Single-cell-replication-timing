load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')
load('data/processed/reference_bulk_profiles.mat', 'ref')
load('data/intermediate/GM12878-unsorted-library1-replicate1.mat', 'g1_windows', ...
    'filtered_counts_20kb', 'is_masked_window')

% Distribution of observed read depths
filtered_counts_20kb = cell2mat(filtered_counts_20kb);
mean_read_depth = nanmean(filtered_counts_20kb);
mean_read_depth = mean_read_depth(~isnan(mean_read_depth));
clearvars filtered_counts_20kb

% Define true states
rt_bins = linspace(-2.3, 2.3, 1000);

bulk = cell(22, 1);
for Chr = 1:22
    bulk{Chr} = interp1(ref.GM12878{Chr}(:, 1), ref.GM12878{Chr}(:, 2), g1_windows{Chr}(:, 3));
end
clearvars ref

true_state = cell(22, 1);
for Chr = 1:22
    true_state{Chr} = double(bulk{Chr} >= rt_bins);
    true_state{Chr}(isnan(bulk{Chr}), :) = NaN;
    true_state{Chr}(true_state{Chr} == 1) = 4;
    true_state{Chr}(true_state{Chr} == 0) = 2;
    true_state{Chr} = true_state{Chr}(:, end:-1:1);
end

flat_states = cell2mat(true_state);
percent_replicated = sum(flat_states == 4, 1) ./ sum(~isnan(flat_states), 1);

% Interpolate true states into fixed windows
true_20kb = cell(22, 1);
for Chr = 1:22
    true_20kb{Chr} = interp1(g1_windows{Chr}(:, 3), true_state{Chr}, genome_windows{Chr}(:, 3));
    is_valid_cn = ismember(true_20kb{Chr}, [2 4]);
    true_20kb{Chr}(~is_valid_cn) = NaN;
end
clearvars bulk Chr rt_bins is_valid_cn flat_states

% Sample with replacement
num_cells = 2500;

rng(50)
index_sim = randsample(1:size(true_state{1}, 2), num_cells, true);
percent_replicated = percent_replicated(index_sim);

rng(23)
mean_read_depth = randsample(mean_read_depth, num_cells, true);

for Chr = 1:22
    true_state{Chr} = true_state{Chr}(:, index_sim);
    true_20kb{Chr} = true_20kb{Chr}(:, index_sim);
end

seeds = cell(22, 1);
for Chr = 1:22
    rng(Chr+42)
    seeds{Chr} = randi(1e6, [num_cells 2]);
end

% Simulate raw counts
raw_counts_20kb = cell(22, 1);
for Chr = 1:22
    disp(['Chromosome ' num2str(Chr)])
    
    raw_counts_20kb{Chr} = NaN(size(true_state{Chr}));
    for barcode = 1:num_cells
                
        cn2 = true_state{Chr}(:, barcode) == 2;
        cn4 = true_state{Chr}(:, barcode) == 4;
        
        
        cn = mean_read_depth(barcode) / ...
            (2* percent_replicated(barcode) + 1 - percent_replicated(barcode));
        
        rng(seeds{Chr}(barcode, 1))
        raw_counts_20kb{Chr}(cn2, barcode) = poissrnd(cn, [sum(cn2) 1]);
        
        rng(seeds{Chr}(barcode, 2))
        raw_counts_20kb{Chr}(cn4, barcode) = poissrnd(2 * cn, [sum(cn4) 1]);
    end
end
clearvars barcode Chr cn cn2 cn4 seeds index_sim

%% Run HMM

% Keep track of chromosome boundaries
chrom_boundaries = NaN(22, 1);
for Chr = 1:22
    chrom_boundaries(Chr) = size(g1_windows{Chr}, 1);
end
chrom_boundaries = cumsum(chrom_boundaries);
chrom_boundaries = [1; chrom_boundaries(1:end-1) + 1];

% Infer replication state by HMM
flat_raw_counts = cell2mat(raw_counts_20kb);
num_cells = size(flat_raw_counts, 2);
flat_masked_windows = cell2mat(is_masked_window);

flat_replication_inference = NaN(size(flat_raw_counts));

parfor barcode = 1:num_cells
    disp(['Cell ' num2str(barcode)])
    [~, flat_replication_inference(:, barcode)] = ...
        infer_single_cell_states(flat_raw_counts(:, barcode), chrom_boundaries);   
end

flat_replication_inference(flat_masked_windows, :) = NaN;

% Un-flatten chromosomes
chrom_boundaries(:, 2) = [chrom_boundaries(2:end, 1) - 1; size(flat_replication_inference, 1)];

replication_inference_g1_windows = cell(22, 1);
for Chr = 1:22
    index = chrom_boundaries(Chr, :);
    replication_inference_g1_windows{Chr} = flat_replication_inference(index(1):index(2), :);
end

% Interpolate to fixed-sized windows
inferred_state = cell(22, 1);
for Chr = 1:22
    inferred_state{Chr} = interp1(g1_windows{Chr}(~is_masked_window{Chr}, 3), ...
        replication_inference_g1_windows{Chr}(~is_masked_window{Chr}, :), ...
        genome_windows{Chr}(:, 3));
    inferred_state{Chr}(~is_mappable{Chr}, :) = NaN;
    inferred_state{Chr}(~ismember(inferred_state{Chr}, 2:4)) = NaN;
end

% Fill across gaps if the gap is smaller than 100Kb AND the state is the same on both sides

for Chr = 1:22
    
    m = find(~is_mappable{Chr});
    index = [1; find(diff(m) > 1) + 1];
    index(:, 2) = [index(2:end) - 1; size(m, 1)];
    is_missing = [m(index(:, 1)) m(index(:, 2))];
    is_missing(:, 3) = is_missing(:, 2) - is_missing(:, 1) + 1;
    is_missing = is_missing(is_missing(:, 3) <= 5, 1:2);
    
    if isempty(is_missing)
        break
    end
    
    for row = 1:size(is_missing, 1)
        range = is_missing(row, 1)-1:is_missing(row, 2)+1;
        
        if any(range <= 0 | range >= size(is_mappable{Chr}, 1))
            continue
        end
        
        values = inferred_state{Chr}(range, :);
        to_fill = values(1, :) == values(end, :);
        for inner = 2:size(values, 1) - 1
            values(inner, to_fill) = values(1, to_fill);
        end
        
        inferred_state{Chr}(range, :) = values;
    end
    
end

clearvars Chr chrom_boundaries flat_masked_windows flat_raw_counts flat_replication_inference ...
    index inner is_missing m range row replication_inference_g1_windows to_fill values ...
    g1_windows is_mappable is_masked_window raw_counts_20kb true_state

%% Assess the impact of 4N region size on calling

% Filter out cells that are >90% incorrect (i.e. G1 cells called as G2 cells)
flat_inference = cell2mat(inferred_state);
p = sum(flat_inference == 4, 1) ./ sum(~isnan(flat_inference), 1);
is_bad_cell = find(abs(p - percent_replicated) > 0.9 | percent_replicated < 0.03);

% Identify missed origins
n_missed = NaN(100, 22);
n_segments = NaN(100, 22);
n_overall = NaN(2, 22);

parfor Chr = 1:22
    
    disp(['Chromosome ' num2str(Chr)])
    
    true_segments = list_cn_segments(true_20kb{Chr}, 4);
    true_segments(:, 4) = true_segments(:, 3) - true_segments(:, 2) + 1;
    
    is_missed = false(size(true_segments, 1), 1);
    for row = 1:size(true_segments, 1)
        barcode = true_segments(row, 1);
        cn = inferred_state{Chr}(true_segments(row, 2):true_segments(row, 3), barcode);
        is_missed(row) = any(cn == 2) & ~any(cn == 4);
        true_segments(row, 5) = sum(isnan(cn)) / length(cn);
    end
    
    is_missed(ismember(true_segments(:, 1), is_bad_cell)) = false;
    is_missed(true_segments(:, 5) > 0.25) = false;
    
    for segment_length = 1:100
        is_length = true_segments(:, 4) == segment_length;
        n_missed(segment_length, Chr) = sum(is_missed(is_length));
        n_segments(segment_length, Chr) = sum(is_length);
    end
    
    n_overall(:, Chr) = [sum(is_missed) size(is_missed, 1)];
    
end

n_missed = sum(n_missed, 2);
n_segments = sum(n_segments, 2);
p_missed_by_size = n_missed ./ n_segments .* 100;

n_overall = sum(n_overall, 2);
p_missed_overall = n_overall(1)/n_overall(2) * 100;

%% Save variables

save('-v7.3', 'data/processed/simulated_GM12878.mat', 'true_20kb', 'inferred_state', ...
    'percent_replicated', 'mean_read_depth', 'num_cells', 'p_missed_by_size', 'p_missed_overall')
