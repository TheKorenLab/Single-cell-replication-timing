function mu = interpolated_median(counts)

    naive_median = nanmedian(counts, 1);
    relative_to_median = sign(counts - naive_median) + 1;

    mu = NaN(1, size(relative_to_median, 2));
    for column = 1:size(relative_to_median, 2)
        
        num_counts = histcounts(relative_to_median(:, column), 0:3);
        
        if length(num_counts) ~= 3 || num_counts(2) == 0
            offset = 0;
        else
            offset = ((num_counts(3) - num_counts(1)) / num_counts(2)) / 2;
        end
        
        mu(column) = naive_median(column) + offset;
    end

end
