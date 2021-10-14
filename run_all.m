clear; clc
load_single_cell_project

%% Processing

% Define 20kb windows
run('processing/define_fixed_windows.m')
run('processing/define_fixed_window_mappability.m')

% Process each sequencing library separately
if ~isfolder('data/intermediate')
    mkdir('data/intermediate')
end

file_list = dir('data/raw/*.h5');

for l = 1:size(file_list, 1)
    
    sample_name = strsplit(file_list(l).name, '.');
    sample_name = sample_name{1};
    disp(sample_name)
    
    count_reads_in_G1_windows([file_list(l).folder '/' file_list(l).name], sample_name)
    assign_replication_states(sample_name)

end

% Combine replicating cells across replicates for the same cell line
run('processing/combine_sc_replicates.m')

% Calculate S/G1 aggregate profile (~2.5 hours)
file_list = dir('data/intermediate/*.mat');
for f = 1:size(file_list, 1)
    sample_name = strsplit(file_list(f).name, '-');
    file_list(f).sample = sample_name{1};
    prefix = strsplit(file_list(f).name, '.');
    file_list(f).prefix = prefix{1};
end

for sample = 1:length(samples)
    disp(samples{sample})
    input_files = {file_list(strcmp({file_list.sample}, samples{sample})).prefix};
    aggregate_S_G1 = calculate_S_G1_aggregate(input_files);
    save('-append', ['data/processed/' samples{sample} '.mat'], 'aggregate_S_G1')
end
clearvars file_list f sample_name prefix sample input_files aggregate_S_G1

% Aggregate vs. bulk
run('processing/save_bulk_profiles.m')
run('processing/calculate_r_bulk.m')

% Partition GM12878 cells into sub-S phase aggregates
run('processing/partition_GM12878_into_subS_fractions.m')

%% Analysis

load('data/hg37_genome_metadata.mat', 'genome_windows')

for sample = 1:10
    
    disp(samples{sample})
    
    filename = ['data/processed/' samples{sample} '.mat'];
    load(filename, 'aggregate_S_G1', 'replication_state_filtered', 'percent_replicated_filtered')
    
    % Call high confidence peaks in the aggregate profile
    high_confidence_peaks = call_high_confidence_peaks(aggregate_S_G1);

    % Identify replication tracks (~11 min)
    replication_tracks = identify_replication_tracks(replication_state_filtered);

    % First-approximation of origin locations
    window_IR_frequency = cell(22, 1);
    for Chr = 1:22
        is_included = replication_tracks{Chr}(:, 4) == 0 & replication_tracks{Chr}(:, 6) <= 50;
        tracks = replication_tracks{Chr}(is_included, 5);
        window_IR_frequency{Chr}(:, 1) = histcounts(tracks, 0:size(genome_windows{Chr}, 1));
    end

    % Define initiation regions (~12 min)
    [replication_tracks, single_cell_IRs] = call_initiation_regions(replication_tracks, ...
        aggregate_S_G1, high_confidence_peaks);

    % Assign tracks to IRs
    [replication_tracks, single_cell_IRs] = assign_tracks_to_IRs(replication_tracks, ...
        single_cell_IRs);

    % Analyze IR firing order
    firing_order_frequency = calculate_out_of_order_firing(single_cell_IRs, ...
        replication_state_filtered);

    % Analyze IR firing time constraint
    [IR_range, cumulative_earliest_timing] = calculate_IR_range(single_cell_IRs, ...
        replication_state_filtered, replication_tracks, percent_replicated_filtered);

    % Classify IRs
    [p_cells_fired, late_IR_class] = classify_late_IRs(single_cell_IRs, ...
        replication_state_filtered, percent_replicated_filtered, replication_tracks);

    % Mask replication states
    replication_state_masked = mask_replication_state(replication_state_filtered, ...
        replication_tracks);

    save('-append', filename, 'window_IR_frequency', 'replication_tracks', ...
        'single_cell_IRs', 'firing_order_frequency', 'IR_range', 'cumulative_earliest_timing', ...
        'p_cells_fired', 'late_IR_class', 'replication_state_masked')

end

% Aggregate profile from G1 fraction
aggregate_G1_fraction = calculate_S_G1_aggregate({'GM12878-G1-library1-replicate1'});
save('-append', 'data/processed/GM12878.mat', 'aggregate_G1_fraction')

% Cell line PCA
run('analysis/run_cell_line_pca.mat')

%% Visualization

run('visualization/figure1_in_silico_sorting.m')
run('visualization/figure2_10x_vs_DLP.m')
run('visualization/figure3_calling_IRs.m')
run('visualization/figure4_IR_variability.m')
run('visualization/figure5_late_IR_behaviors.m')
run('visualization/figure6_additional_cell_lines.m')
run('visualization/figureS1_mapd_validation.m')
run('visualization/figureS2_wide_IRs.m')
run('visualization/figureS3_cell_line_QC.m')
run('visualization/figureS4_chromosome_plots.m')
run('visualization/figureS5_cell_line_specificity.m')
run('visualization/figureS6_7_cell_line_variability.m')
