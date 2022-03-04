load_single_cell_project
load('data/processed/GM12878.mat', 'aggregate_S_G1', 'single_cell_IRs')

aggregate_S_G1 = cell2mat(aggregate_S_G1);
aggregate_S_G1 = aggregate_S_G1(:, 2);

single_cell_IRs = cell2mat(single_cell_IRs);
single_cell_IRs = single_cell_IRs(:, 11);

[rt, edges] = histcounts(aggregate_S_G1, 20);
rt = rt ./ sum(rt) .* 100;

[irs, ~] = histcounts(single_cell_IRs, edges);
irs = irs ./ sum(irs) .* 100;

x = [edges(1:end-1); edges(2:end)];
x = mean(x, 1);

figureR5 = figure;
set(figureR5, 'Units', 'inches', 'Position', [25 9 3 2.5])
set(gca, 'Position', [0.8848 0.7108 6.4530 5.1629], 'XLim', [-1.9 2.1])
bar(x, rt, 1, 'FaceAlpha', 0.25, 'DisplayName', 'Genome windows')
bar(x, irs, 1, 'FaceAlpha', 0.25, 'DisplayName', 'Initiation regions')
xlabel('Repliction timing')
ylabel('% of total')
title('Replication timing distribution')
lh = legend('show');
lh.ItemTokenSize(1) = 10;

printFigure('out/FigureR5.pdf')
close all
clearvars
