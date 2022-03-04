load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows', 'is_mappable')

file_list = dir('data/intermediate/*.mat');
for f = 1:size(file_list, 1)
    sample_name = strsplit(file_list(f).name, '-');
    file_list(f).sample = sample_name{1};
    prefix = strsplit(file_list(f).name, '.');
    file_list(f).prefix = prefix{1};
end

% Load single cell data
data = struct;
for sample = 1:10
    
    input_files = {file_list(strcmp({file_list.sample}, samples{sample})).prefix};
    
    mean_coverage_1Mb = cell(1, length(input_files));
    scaled_mapd_1Mb = cell(1, length(input_files));
    is_G1 = cell(1, length(input_files));
    for f = 1:length(input_files)
        filename = ['data/intermediate/' input_files{f} '.mat'];
        in = load(filename, 'mean_coverage_1Mb', 'scaled_mapd_1Mb', 'is_G1');
        mean_coverage_1Mb{f} = in.mean_coverage_1Mb;
        scaled_mapd_1Mb{f} = in.scaled_mapd_1Mb;
        is_G1{f} = in.is_G1;
    end

    data.(samples{sample}).mean_coverage_1Mb = cell2mat(mean_coverage_1Mb);
    data.(samples{sample}).scaled_mapd_1Mb = cell2mat(scaled_mapd_1Mb);
    data.(samples{sample}).is_G1 = cell2mat(is_G1);
    
    clearvars is_G1 mean_coverage_1Mb scaled_mapd_1Mb in
end

% Load bulk seq
load('data/processed/reference_bulk_profiles.mat', 'ref')

% MCF-7 fixed vs. G1 windows
hdf_file = 'data/raw/MCF7-unsorted-library1.h5';
coverage_20kb = cell(22, 1);
for Chr = 1:22
    coverage_20kb{Chr} = double(h5read(hdf_file, ['/raw_counts/chr' num2str(Chr)]));
    coverage_20kb{Chr}(~is_mappable{Chr}, :) = NaN;
end

load('data/intermediate/MCF7-unsorted-library1.mat', 'raw_counts_20kb', 'g1_windows')

data.MCF7.g1_windows = g1_windows;
data.MCF7.fixed_window_counts = coverage_20kb;
data.MCF7.g1_window_counts = raw_counts_20kb;

clearvars cancer_bulk Chr coverage_20kb ESC_consensus hdf_file index is_mappable LCL_bulk ...
    raw_counts_20kb sample

%% Figure skeleton

figureS6 = figure;
set(figureS6, 'Units', 'inches', 'Position', [25 9 6.5 5.85])

x = linspace(0.45, 5.4, 5);
y = [4.51 2.84];

s = 0;
for yi = 1:2
    for xi = 1:5
        s = s + 1;
       
        if s > 10
            continue
        end
        
        panelA.(samples{s}) = axes('Units', 'inches', 'Position', [x(xi) y(yi) 1 0.89]);
    end
end

panelB = struct('top', axes('Units', 'inches', 'Position', [0.45 0.99 5.8 0.75]), ...
    'bottom', axes('Units', 'inches', 'Position', [0.45 0.44 5.8 0.5]));

%% Panel A

params = struct('color', {g1_dark s_dark}, 'labels', {'G1', 'S'}, 'alpha', {1 0.5});

for p = 1:10
    
    parent = panelA.(samples{p});
    
    fraction = [data.(samples{p}).is_G1; ~data.(samples{p}).is_G1]';
    
    for f = 1:2
        scatter(data.(samples{p}).mean_coverage_1Mb(fraction(:, f)), ...
            data.(samples{p}).scaled_mapd_1Mb(fraction(:, f)), '.', ...
            'MarkerFaceColor', params(f).color, 'MarkerFaceAlpha', params(f).alpha, ...
            'MarkerEdgeColor', params(f).color, 'MarkerEdgeAlpha', params(f).alpha, ...
            'HandleVisibility', 'off', 'Parent', parent)
    end
    
    Xmax = prctile(data.(samples{p}).mean_coverage_1Mb, 99);
    XTick = round(linspace(0, Xmax-20, 4) / 5) .* 5;
    Ymax = prctile(data.(samples{p}).scaled_mapd_1Mb, 99);
    set(parent, 'XLim', [20 Xmax], 'XTick', XTick, 'YLim', [0.8 Ymax], 'YTick', [0 1 2])
        
    xlabel(parent, 'Reads per Mb')
    
    if ismember(p, [1 6])
        ylabel(parent, 'Scaled MAPD')
    end
    
    title(parent, cell_line_names{p})
end

% Legend
legend_markers = cell(2, 1);
for f = 1:2
    legend_markers{f} = scatter(0, 0, 20, 'filled', 'Parent', panelA.GM12878);
end
set(legend_markers{1}, 'MarkerFaceColor', g1_dark, 'DisplayName', 'G1')
set(legend_markers{2}, 'MarkerFaceColor', s_dark, 'DisplayName', 'S')

legendA = legend(panelA.GM12878);
legendA.ItemTokenSize(1) = 5;
set(legendA, 'Orientation', 'horizontal', 'Units', 'inches', 'Position', [3.08 5.6148 0.65 0.1667])

%% Panel B

params = struct('Chr', 1, 'panel', {panelB.top panelB.bottom}, ...
    'data', {data.MCF7.fixed_window_counts data.MCF7.g1_window_counts}, ...
    'windows', {genome_windows data.MCF7.g1_windows}, 'label', {'fixed', 'G1'});

for p = 1:2
    
    parent = params(p).panel;
    
    d = sum(params(p).data{params(p).Chr}(:, ~data.MCF7.is_G1), 2);
    d = 2 .* d ./ nanmedian(d);
    
    plot(params(p).windows{params(p).Chr}(:, 3) ./ 1e6, d, '.', 'Color', s_light{3}, ...
        'MarkerSize', 3, 'Parent', parent)
    
    X = [genome_windows{params(p).Chr}(1, 1) genome_windows{params(p).Chr}(end, 2)]./1e6;
    set(parent, 'XLim', X)
    ylabel(parent, 'CN');
    
    if p == 1
        set(parent, 'XTick', [], 'YTick', [2 10 20])
    else
        set(parent, 'YLim', [0.75 3.25])
    end
    
    yh = get(parent,'YLabel');
    yh.Position(1) = -10;
    
    yyaxis(parent, 'right')
    set(parent, 'YColor', 'k', 'YTick', [])
    ylabel(parent, params(p).label)
    
end
    
title(panelB.top, 'MCF-7, fixed-size windows vs. G1-windows')
xlabel(panelB.bottom, ['Chromosome ' num2str(params(p).Chr) ' Coordinate, Mb'])

%% Annotate panels

params = struct('panel', {panelA.GM12878 panelB.top}, 'text', {'a', 'b'}, 'X', {-0.2361 -0.0479}, ...
    'Y', {1.2171 1.3333});

for p = 1:2
    text(params(p).X, params(p).Y, params(p).text, 'Parent', params(p).panel, ...
        'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS6.pdf')
close
clearvars
