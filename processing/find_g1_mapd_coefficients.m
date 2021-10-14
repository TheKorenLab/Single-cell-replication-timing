function g1_coefficients = find_g1_mapd_coefficients(mean_coverage, scaled_mapd)

    num_cells = length(mean_coverage);
    assert(length(scaled_mapd) == num_cells)

    % Initialize group assignments randomly
    group = randi(2, [1 num_cells]);
    previous_assignments = ones(1, num_cells);

    % Fit two-component linear model
    coefficients = zeros(2, 2);
    residuals = zeros(2, num_cells);
    likelihoods = zeros(2, num_cells);

    while sum(previous_assignments ~= group) > 0
        
        for group_i = 1:2
            coefficients(group_i, :) = polyfit(mean_coverage(group == group_i), ...
                scaled_mapd(group == group_i), 1);
            residuals(group_i, :) = scaled_mapd - polyval(coefficients(group_i, :), mean_coverage);
            likelihoods(group_i, :)  = normpdf(residuals(group_i, :), 0, std(residuals(group_i, :)));
        end

        previous_assignments = group;
        for barcode = 1:num_cells
            group(barcode) = find(likelihoods(:, barcode) == max(likelihoods(:, barcode)));
        end

    end

    g1_coefficients = coefficients(coefficients(:, 2) == min(coefficients(:, 2)), :);

    if isempty(g1_coefficients)
        error('G1 population was not detected successfully.')
    end

end
