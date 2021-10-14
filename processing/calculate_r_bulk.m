load_single_cell_project
load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/reference_bulk_profiles.mat', 'ref')

data = struct;
for sample = 1:10
    data.(samples{sample}) = load(['data/processed/' samples{sample} '.mat'], 'aggregate_S_G1');
end

r_bulk = struct;
for sample = 1:10

    profiles = cell(22, 1);
    for Chr = 1:22
        profiles{Chr}(:, 1) = interp1(ref.(samples{sample}){Chr}(:, 1), ...
            ref.(samples{sample}){Chr}(:, 2), genome_windows{Chr}(:, 3));
        
        index = ~isnan(data.(samples{sample}).aggregate_S_G1{Chr}(:, 2));
        profiles{Chr}(:, 2) = interp1(data.(samples{sample}).aggregate_S_G1{Chr}(index, 1), ...
           data.(samples{sample}).aggregate_S_G1{Chr}(index, 2), genome_windows{Chr}(:, 3)); 
    end
    
    profiles = cell2mat(profiles);
    r_bulk.(samples{sample}) = corr(profiles(:, 1), profiles(:, 2), 'rows', 'pairwise');
end

save('-append', 'data/processed/reference_bulk_profiles.mat', 'r_bulk')
