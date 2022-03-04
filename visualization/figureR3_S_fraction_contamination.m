load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')

hdf_file = 'data/raw/GM12878-S-library1-replicate1.h5';
coverage_20kb = cell(22, 1);
for Chr = 1:22
    coverage_20kb{Chr} = double(h5read(hdf_file, ['/raw_counts/chr' num2str(Chr)]));
    coverage_20kb{Chr}(~is_mappable{Chr}, :) = NaN;
end

load('data/intermediate/GM12878-S-library1-replicate1.mat', 'is_G1')

profiles = cell(22, 1);
for Chr = 1:22
    profiles{Chr}(:, 1) = nansum(coverage_20kb{Chr}(:, is_G1), 2);
    profiles{Chr}(:, 2) = nansum(coverage_20kb{Chr}(:, ~is_G1), 2);
    profiles{Chr}(profiles{Chr} == 0) = NaN;
end

Chr = 1;
X = [0 120];

figureR3 = figure;
set(figureR3, 'Units', 'inches', 'Position', [25 12 6.5 2.67])

params = struct('y', {1.65 0.35}, 'ymax', {1500 5000}, 'title', {'G1/G2', 'S'});

for p = 1:2
    parent = axes('Units', 'inches', 'Position', [0.4583 params(p).y 5.89 0.8]);
    
    plot(genome_windows{Chr}(:, 3) ./ 1e6, profiles{Chr}(:, p), 'k.', 'MarkerSize', 3, ...
        'Parent', parent)
    set(parent, 'XLim', X, 'YLim', [0 params(p).ymax])
    xlabel(parent, ['Chromosome ' num2str(Chr) ' Coordinate, Mb'])
    ylabel(parent, 'Count')
    title(parent, [params(p).title '-phase cells in the S-phase Fraction'])

end

printFigure('out/FigureR3.pdf')
close
clearvars

