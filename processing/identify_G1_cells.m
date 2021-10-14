function [mean_coverage, scaled_mapd, is_G1, g1_barcodes] = identify_G1_cells(coverage_20kb, ...
    is_mappable)

    % Identify G1 cells by MAPD
    coverage_1Mb = cell(22, 1);
    for Chr = 1:22
        counts = coverage_20kb{Chr}(is_mappable{Chr}, :);
        coverage_1Mb{Chr} = aggregate_counts(counts, 50);
    end

    flat_coverage_1Mb = cell2mat(coverage_1Mb);

    mean_coverage = nanmean(flat_coverage_1Mb, 1);
    is_low_coverage = mean_coverage <= 10;

    scaled_mapd = calculate_scaled_mapd(flat_coverage_1Mb, mean_coverage);
    g1_coefficients = find_g1_mapd_coefficients(mean_coverage(~is_low_coverage), ...
        scaled_mapd(~is_low_coverage));

    residual_mapd = scaled_mapd - polyval(g1_coefficients, mean_coverage);
    is_G1 = residual_mapd < 0.05;
    is_G1(is_low_coverage) = false;

    % Identify high confidence G1 cells to use for windowing
    % Sequentially include more cells until 200 reads ~ 20kb

    [~, barcodes] = sort(mean_coverage, 'descend');
    barcodes = barcodes(is_G1(barcodes));

    flat_coverage_20kb = cell2mat(coverage_20kb);
    flat_coverage_20kb(flat_coverage_20kb == 0) = NaN;
 
    num_G1_cells = 0;
    window_size = 200 / nanmean(nansum(flat_coverage_20kb(:, barcodes(1)), 2), 1);
    while window_size > 1 && num_G1_cells < sum(is_G1)
        num_G1_cells = num_G1_cells + 1;
        window_size = 200 / ...
            nanmean(nansum(flat_coverage_20kb(:, barcodes(1:num_G1_cells)), 2), 1);
    end

    g1_barcodes = sort(barcodes(1:num_G1_cells), 'ascend');

end
