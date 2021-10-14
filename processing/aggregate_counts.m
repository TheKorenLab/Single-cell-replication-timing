function aggregated = aggregate_counts(counts, windows_to_aggregate)

    aggregated_height = floor((size(counts, 1)-1)/windows_to_aggregate) + 1;
    
    aggregated = zeros(aggregated_height, size(counts, 2));
    index = 1:windows_to_aggregate:size(counts, 1);
    for row = 1:length(index)
        high = min(index(row) + windows_to_aggregate, size(counts, 1) + 1);
        aggregated(row, :) = nansum(counts(index(row):high-1, :), 1);
    end
end
