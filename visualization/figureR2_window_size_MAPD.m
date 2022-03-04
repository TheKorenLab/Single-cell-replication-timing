load_single_cell_project
load('data/hg37_genome_metadata.mat', 'is_mappable')

samples = {'G1', 'earlyS', 'S', 'lateS', 'G2'};
widths = [10 25 50]; % x 20kb

num_G1 = NaN(5, 3);
num_cells = NaN(5, 1);

for s = 1:length(samples)
    
    disp(samples{s})
    
    hdf_file = ['data/raw/GM12878-' samples{s} '-library1-replicate1.h5'];
    coverage_20kb = cell(22, 1);
    for Chr = 1:22
        coverage_20kb{Chr} = double(h5read(hdf_file, ['/raw_counts/chr' num2str(Chr)]));
        coverage_20kb{Chr}(~is_mappable{Chr}, :) = NaN;
    end
   
    num_cells(s) = size(coverage_20kb{1}, 2);

    % Identify G1 cells by MAPD
    is_G1 = false(num_cells(s), 3);
    for w = 1:3
        coverage_aggregated = cell(22, 1);
        for Chr = 1:22
            counts = coverage_20kb{Chr}(is_mappable{Chr}, :);
            coverage_aggregated{Chr} = aggregate_counts(counts, widths(w));
        end

        flat_coverage = cell2mat(coverage_aggregated);
        mean_coverage = nanmean(flat_coverage, 1);
        scaled_mapd = calculate_scaled_mapd(flat_coverage, mean_coverage);

        is_low_coverage = mean_coverage <= 10;
        g1_coefficients = find_g1_mapd_coefficients(mean_coverage(~is_low_coverage), ...
            scaled_mapd(~is_low_coverage));
        residual_mapd = scaled_mapd - polyval(g1_coefficients, mean_coverage);
        
        is_G1(:, w) = residual_mapd < 0.05;
        is_G1(is_low_coverage, w) = false;
    end
    
    num_G1(s, :) = sum(is_G1, 1);
end

num_G1 = num_G1 ./ num_cells .* 100;

figureR2 = figure;
set(figureR2, 'Units', 'inches', 'Position', [25 12 3 2])
set(gca, 'Position', [0.9906 0.7123 6.3472 3.9867])

bar(num_G1, 'grouped')
set(gca, 'XLim', [0.5 5.5], 'XTick', 1:5, 'YLim', [0 110], ...
    'XTickLabel', {'G1', 'Early S', 'S', 'Late S', 'G2'})
xlabel('FACS Fraction')
ylabel('% G1 Cells')
title('Effect of window size')
lh = legend('show', {'200Kb', '500Kb', '1Mb'}, 'Orientation', 'Horizontal');
lh.ItemTokenSize(1) = 15;

printFigure('out/FigureR2.pdf')
close
clearvars
