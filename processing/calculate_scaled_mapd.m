function scaled_mapd = calculate_scaled_mapd(read_counts, mean_coverage)

    absolute_difference = abs(diff(read_counts) - nanmedian(diff(read_counts)));
    mapd = interpolated_median(absolute_difference) ./ mean_coverage;
    scaled_mapd = mapd .* sqrt(mean_coverage);
    
end
