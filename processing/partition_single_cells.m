function subS_fractions = partition_single_cells(counts, n_fractions)

    n_cells = size(counts{1}, 2);

    subS_fractions = cell(22, length(n_fractions));

    for f = 1:length(n_fractions)

        % Aggregate counts
        for Chr = 1:22
            n = floor(n_cells / n_fractions(f)) * n_fractions(f);
            x = floor(linspace(1, n, n_fractions(f) + 1));
            x = [x(1:end-1); x(2:end)]';
            x(2:end, 1) = x(2:end, 1) + 1;

            subS_fractions{Chr, f} = NaN(size(counts{Chr}, 1), n_fractions(f));

            for col = 1:n_fractions(f)
                subS_fractions{Chr, f}(:, col) = nansum(counts{Chr}(:, x(col,1):x(col,2)), 2);
            end
        end

        % Normalize to min 2 and max 4
        flat_counts = cell2mat(subS_fractions(:, f));
        flat_counts(flat_counts == 0) = NaN;
        m = [min(flat_counts, [], 1); prctile(flat_counts, 99.9)];

        for Chr = 1:22
            subS_fractions{Chr, f} = 2 + ((subS_fractions{Chr, f} - m(1, :)) * (4 - 2)) ./ ...
                (m(2, :) - m(1, :));

            subS_fractions{Chr, f}(subS_fractions{Chr, f} < 2) = 2;
            subS_fractions{Chr, f}(subS_fractions{Chr, f} > 4) = 4;
        end

    end

end
