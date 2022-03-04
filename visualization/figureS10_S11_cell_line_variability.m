load_single_cell_project

data = struct;
all_widths = cell(9, 1);
for sample = 1:9
    
    filename = ['data/processed/' samples{sample} '.mat'];
    load(filename, 'single_cell_IRs', 'firing_order_frequency', 'IR_range', ...
        'cumulative_earliest_timing', 'p_cells_fired', 'late_IR_class')
    
    flat_IR_list = cell2mat(single_cell_IRs);
    
    % IR width distribution
    IR_widths = flat_IR_list(flat_IR_list(:, 10) >= 5, 7);
    all_widths{sample} = IR_widths;
    mu_width = nanmedian(IR_widths);
    
    % IR order analysis
    flat_sorted_rt_values = sort(flat_IR_list(:, 11), 'descend');
    firing_order_frequency(:, 5) = sum(firing_order_frequency(:, 3:4), 2);
    firing_order_frequency = [flat_sorted_rt_values firing_order_frequency]; %#ok<AGROW>
    
    index = ~isnan(flat_sorted_rt_values);
    XTicks = flip(interp1(flat_sorted_rt_values(index), 1:sum(index), [-1.25 0 1.25]));
    
    % IR firing time analysis
    IR_range = [flat_IR_list(:, 11) cell2mat(IR_range)];
    IR_range(:, 2:3) = IR_range(:, 2:3) .* 100;
    IR_range(any(isnan(IR_range(:, 2:3)), 2), 4) = 0;
    IR_range(:, 5) = IR_range(:, 3) - IR_range(:, 2);
    IR_range = IR_range(:, [1:3 5 4]);
    IR_range = sortrows(IR_range, 1, 'descend');
    
    % Late IR classes analysis
    p_cells_fired = cell2mat(p_cells_fired);
    late_IR_class = cell2mat(late_IR_class);
    late_IR_class = late_IR_class(:, [3 2 1]);
    
    data.(samples{sample}) = struct('IR_widths', IR_widths, 'mu_width', mu_width, ...
        'firing_order_frequency', firing_order_frequency, 'IR_range', IR_range, 'XTicks', XTicks, ...
        'cumulative_earliest_timing', cumulative_earliest_timing, ...
        'p_cells_fired', p_cells_fired, 'late_IR_class', late_IR_class);
    
end

% Consistent edges for histograms
all_widths = cell2mat(all_widths);
[~, IR_width_bins] = histcounts(all_widths / 1e3, 25);
IR_width_x = [IR_width_bins(1:end-1); IR_width_bins(2:end)];
IR_width_x = mean(IR_width_x, 1);

clearvars -except data samples IR_width_bins IR_width_x

sample_names = samples(1:9);
sample_names{7} = 'HCT-116';
sample_names{9} = 'MCF-7';

bar_colors = {'#DE2C26', '#3182BD', '#636362', '#BDBDBD'};
range_colors = interp1([0 100], [215 48 39; 145 191 219]./255, linspace(0, 100, 128), 'pchip');
pie_colors = {'#66C2A5'; '#E78AC3'; '#8DA0CB'; '#FC8D62'};

%% Figure S6

figureS6 = figure;
set(figureS6, 'Position', [25 9 6.5 9])

% Figure skeleton
panel_names = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'};

panels = struct;
params = struct('panel', {1 2 3}, 'y', {8.01 4.915 1.8}, ...
    'height', {0.64 0.9 0.66}, 'width', {0.85 1.12 1.1});

for p = 1:3
    
    x = linspace(0.45, 6.4-params(p).width, 5);
    y = [params(p).y params(p).y-params(p).height-0.7];
    
    counter = 0;
    for yi = 1:length(y)
        for xi = 1:length(x)
            counter = counter + 1;
            
            if counter > 9
                continue
            end
            
            panels(p).(panel_names{counter}) = axes('Units', 'inches', ...
                'Position', [x(xi) y(yi) params(p).width params(p).height]);
        end
    end
    
end

% Plot data

for sample = 1:9
    
    s = data.(samples{sample});
    
    % IR width
    parent = panels(1).(panel_names{sample});
    
    counts = histcounts(s.IR_widths / 1e3, IR_width_bins);
    counts = counts ./ sum(counts) .* 100;
    
    bar(IR_width_x, counts, 1, 'FaceColor', '#BDBDBD', 'EdgeColor', '#636363', 'Parent', parent);
    axis(parent, 'tight')
    set(parent, 'XLim', [20 500], 'YLim', [0 55])
    xlabel(parent, 'IR width, kb')
    title(parent, sample_names{sample})
    
    text(325, 44.5, 'Median:', 'Parent', parent, 'FontSize', 9, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'center');
    text(325, 30.75, [num2str(s.mu_width/1e3, '%0.2f') 'kb'], ...
        'Parent', parent, 'FontSize', 9, 'FontName', 'Arial', 'HorizontalAlignment', 'center');
    
    if ismember(sample, [1 6])
        ylabel(parent, '% of IRs')
    end
    
    % Predictability of firing order
    parent = panels(2).(panel_names{sample});
    
    b = bar(flip(s.firing_order_frequency(:, 2:5), 2), 1, 'stacked', 'EdgeColor', 'none', ...
        'Parent', parent);
    for row = 1:4
        b(row).FaceColor = bar_colors{row};
    end
    set(parent, 'XLim', [1 size(s.firing_order_frequency, 1)], 'YLim', [0 100], ...
        'YTick', [0 25 75], 'XTick', s.XTicks, 'XTickLabel', [1.25 0 -1.25])
    xlabel(parent, 'Aggregate RT')
    title(parent, sample_names{sample})
    
    if ismember(sample, [1 6])
        ylabel(parent, '% of IRs')
    else
        set(parent, 'YTick', [])
    end
    
    % Out-of-order firing
    
    parent = panels(3).(panel_names{sample});
    
    plot(s.firing_order_frequency(:, 1), s.firing_order_frequency(:, 6), '.', 'MarkerSize', 4, ...
        'Color', '#BDBDBD', 'Parent', parent)
    set(parent, 'XLim', [-2 2], 'XTick', [-1.25 0 1.25], 'YLim', [0 40], 'XDir', 'reverse')
    xlabel(parent, 'Aggregate RT')
    title(parent, sample_names{sample})
    
    if ismember(sample, [1 6])
        ylabel(parent, '% of IRs')
    else
        set(parent, 'YTick', [])
    end
    
end

legendB = legend(panels(2).a, b([4 3 2 1]), {'On Time', 'Yet to Fire', 'Delayed', 'Premature'});
set(legendB,'FontSize', 9, 'Units', 'inches', 'Position', [5.4 3.435 0.88 0.66])
legendB.ItemTokenSize(1) = 10;

% Annotate panels
params = struct('text', {'a', 'b', 'c'}, 'x', {-0.3607 -0.2862 -0.3101}, ...
    'y', {1.2796 1.2021 1.2796});

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', panels(p).a, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS6.pdf')
close

clearvars -except range_colors pie_colors data samples sample_names

%% Figure S7

figureS7 = figure;
set(figureS7, 'Position', [25 9 6.5 9])

% Figure skeleton
panel_names = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'};

panels = struct;
params = struct('panel', {1 2 3}, 'x', {[0.45 4.48] [0.45 2.3168] [3.54 5.6139]}, ...
    'y', {[7.5306 4.3584] [2.5 0.44] [2.375, 0.0644]}, ...
    'height', {1.1 0.61 0.78}, 'width', {1.92 0.72 0.78});

for p = 1:3
    
    x = linspace(params(p).x(1), params(p).x(2), 3);
    y = linspace(params(p).y(1), params(p).y(2), 3);
    
    counter = 0;
    for yi = 1:length(y)
        for xi = 1:length(x)
            counter = counter + 1;
            
            if counter > 9
                continue
            end
            
            panels(p).(panel_names{counter}) = axes('Units', 'inches', ...
                'Position', [x(xi) y(yi) params(p).width params(p).height]);
        end
    end
    
end

% Plot data

for sample = 1:9
    
    s = data.(samples{sample});
    
    % Firing range
    
    parent = panels(1).(panel_names{sample});
    
    index = logical(s.IR_range(:, 5));
    c = interp1(linspace(0, 100, 128), range_colors, s.IR_range(:, 4));
    
    for o = 1:size(s.IR_range, 1)
        if ~index(o)
            continue
        end
        plot(o * ones(1, 2), s.IR_range(o, 2:3), 'Color', c(o, :), 'LineWidth', 0.6, ...
            'Parent', parent)
    end
    set(parent, 'YDir', 'reverse', 'XLim', [1 size(s.IR_range, 1)] + [-5 5], 'CLim', [0 100], ...
        'XTick', s.XTicks, 'XTickLabel', flip(-1.5:0.5:1), 'YTick', 20:20:90)
    title(parent, sample_names{sample})
    
    if ismember(sample, [1 4 7])
        ylabel(parent, '% Replicated')
    else
        set(parent, 'YTick', [])
    end
    
    if sample >= 7
        xlabel(parent, 'Aggregate RT')
    end
    
    % Cumulative early timing
    
    parent = panels(2).(panel_names{sample});
    
    plot(s.cumulative_earliest_timing(:, 1) * 100, s.cumulative_earliest_timing(:, 2) * 100, ...
        '.', 'MarkerSize', 4, 'Color', 'k', 'Parent', parent)
    set(parent, 'XLim', [5 75], 'XTick', [25 50], 'YLim', [0 105], 'YTick', [25 75])
    title(parent, sample_names{sample})
    
    if ismember(sample, [1 4 7])
        ylabel(parent, '% IRs')
    else
        set(parent, 'YTick', [])
    end
    
    if sample >= 7
        xlabel(parent, 'Earliest Time')
    end
    
    % IR class pie chart
    parent = panels(3).(panel_names{sample});
    
    pie_values = [sum(s.p_cells_fired > 0.5) sum(s.late_IR_class, 1)];
    colors = pie_colors(pie_values ~= 0);
    p = pie(pie_values, ones(1,4), 'Parent', parent);
    
    set(parent, 'XColor', 'none', 'YColor', 'none')
    
    t = findobj(p, 'Type', 'text');
    
    for l = 1:length(t)
        if contains(t(l).String, '<')
            set(t(l), 'Visible', 'off')
        end
        
        n = t(l).String;
        n = str2double(n(1:end-1));
        
        if n < 10
            t(l).String = num2str(n);
        end
    end
    
    set(t, 'FontSize', 7.5, 'FontName', 'Arial', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle')
    t_pos = cell2mat(get(t,'Position'));
    set(t, {'Position'}, num2cell(t_pos * 0.68, 2))
    
    c = findobj(p, 'Type', 'patch');
    set(c, {'FaceColor'}, colors, 'LineWidth', 0.6)
    
    title(parent, sample_names{sample})
    
end

legendC = legend(panels(3).a, {'Early in aggregate', 'Throughout S', 'Early + Rare', ...
    'Constitutively Late'});
set(legendC, 'FontSize', 9, 'Units', 'inches', 'Position', [3.7235 3.4306 2.6319 0.4045], ...
    'NumColumns', 2)
legendC.ItemTokenSize(1) = 10;


% Annotate panels
params = struct('text', {'a', 'b', 'c'}, 'x', {-0.1232 -0.4466 -0.1964}, ...
    'y', {1.1509 1.4773 1.3571});

for p = 1:size(params, 2)
    text(params(p).x, params(p).y, params(p).text, 'Parent', panels(p).a, ...
        'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Units', 'normalized');
end

printFigure('out/FigureS7.pdf')
close
