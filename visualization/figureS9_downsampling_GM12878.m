load_single_cell_project

load('data/processed/cell_number_vs_IRs.mat', 'sample_sizes', 'num_IRs')

cell_line_colors = {'#41AE76','#238B45', '#005824', '#74A9CF', '#2B8CBE', '#045A8D', '#F768A1', ...
    '#C51B8A', '#7A0177'};

num_observed = NaN(9, 2);
for sample = 1:9
    filename = ['data/processed/' samples{sample} '.mat'];
    load(filename, 'single_cell_IRs', 'barcodes')
    
    flat = cell2mat(single_cell_IRs);
    num_observed(sample, :) = [length(barcodes) size(flat, 1)];
end

figureS9 = figure;
set(figureS9, 'Units', 'inches', 'Position', [25 12 4.5 4.5])

for sample = 1:9
    plot(-20, 0, 'o', 'MarkerEdgeColor', cell_line_colors{sample}, ...
        'MarkerFaceColor', cell_line_colors{sample}, 'MarkerSize', 6, ...
        'DisplayName', cell_line_names{sample})
    
    plot(num_observed(sample, 1), num_observed(sample,2),'o', ...
        'MarkerEdgeColor', cell_line_colors{sample}, ...
    'MarkerFaceColor', cell_line_colors{sample}, 'MarkerSize', 6, 'HandleVisibility', 'off')

end

plot(-20, 0, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', ...
    'DisplayName', 'Downsampled GM12878')

for n = 1:length(sample_sizes)
    plot(sample_sizes(n), num_IRs(:, n), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', ...
        'HandleVisibility', 'off')
end


xlabel('Cells sampled')
ylabel('IRs called')
title('Effects of Sample Size on IR Calling')
xlim([0 2800])

lh = legend('show');
set(lh, 'Units', 'inches', 'Position', [2.6980 0.5764 1.2986 1.2292])
lh.ItemTokenSize(1) = 5;

printFigure('out/FigureS9.pdf')
close
clearvars
