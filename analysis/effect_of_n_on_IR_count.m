load_single_cell_project
load('data/hg37_genome_metadata.mat', 'genome_windows')

sample = 1;
filename = ['data/processed/' samples{sample} '.mat'];
load(filename, 'aggregate_S_G1', 'replication_state_filtered')

% Call high confidence peaks in the aggregate profile
high_confidence_peaks = call_high_confidence_peaks(aggregate_S_G1);

% Identify replication tracks
replication_tracks = identify_replication_tracks(replication_state_filtered);

sample_sizes = [100 500 1000 2000];

barcodes = 1:size(replication_state_filtered{1}, 2);

rng(10)
seeds = randi(max(barcodes), [5 4]);

tic
num_IRs = NaN(size(seeds));
for n = 1:length(sample_sizes)
    for it = 1:5
        rng(seeds(it, n))
        sample = sort(randsample(barcodes,  sample_sizes(n), false));
        
        r = cell(22, 1);
        for Chr = 1:22
            index = ismember(replication_tracks{Chr}(:, 1), sample);
            r{Chr} = replication_tracks{Chr}(index, :);
        end

        % Define initiation regions
        [r, single_cell_IRs] = call_initiation_regions(r, aggregate_S_G1, high_confidence_peaks);
        single_cell_IRs = cell2mat(single_cell_IRs);
        num_IRs(it, n) = size(single_cell_IRs, 1);
    end
end
toc

save('data/processed/cell_number_vs_IRs.mat', 'num_IRs', 'sample_sizes')
