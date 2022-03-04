load_single_cell_project

load('data/hg37_genome_metadata.mat', 'genome_windows')
load('data/processed/GM12878.mat', 'single_cell_IRs', 'late_IR_class', 'replication_tracks', ...
    'percent_replicated_filtered', 'replication_state_filtered')

x = 0:0.1:1;

percent_fired = cell(22, 1);
for Chr = 1:22
       
    percent_fired{Chr} = NaN(length(single_cell_IRs{Chr}), 10);
    for o = 1:size(late_IR_class{Chr}, 1)
        
        if ~late_IR_class{Chr}(o, 1)
            continue
        end
        
        
        rep_cells = replication_tracks{Chr}(replication_tracks{Chr}(:, 8) == o, 1);        
        p_rep_cells = percent_replicated_filtered(rep_cells);

        IR_center = find(genome_windows{Chr}(:, 1) < single_cell_IRs{Chr}(o, 4) & ...
            genome_windows{Chr}(:, 2) >= single_cell_IRs{Chr}(o, 4));
        unrep_cells = replication_state_filtered{Chr}(IR_center, :) == 2;
        p_unrep_cells = percent_replicated_filtered(unrep_cells);
        
        p_informative = [p_rep_cells p_unrep_cells];
        
        for xi = 1:length(x)-1
            n_firing =  sum(p_rep_cells > x(xi) & p_rep_cells <= x(xi+1));
            n_informative = sum(p_informative > x(xi) & p_informative <= x(xi+1));
            percent_fired{Chr}(o, xi) = n_firing / n_informative;
        end

    end
end

% IR number
for Chr = 2:22
    single_cell_IRs{Chr}(:, 1) = single_cell_IRs{Chr}(:, 1) + single_cell_IRs{Chr-1}(end, 1);
end
IR_number = cell2mat(single_cell_IRs);
IR_number = IR_number(:, 1);

% Filter out the other clasess of IRs
percent_fired = cell2mat(percent_fired);
index = all(isnan(percent_fired), 2);
percent_fired = percent_fired(~index, :);
IR_number = IR_number(~index);

% Figure
x = linspace(0.38, 5.47, 5);
y = linspace(4.8258, 0.3128, 4);
pos = cell(5, 4);
for xi = 1:5
    for yi = 1:4
        pos{xi, yi} = [x(xi) y(yi)];
    end
end

pos = cell2mat(reshape(pos, [20 1]));

rng(52313)
irs_shown = sort(randsample(1:size(percent_fired, 1), 20));
ymax = max(max(percent_fired(irs_shown, :))) * 100;

figureR6 = figure;
set(figureR6, 'Units', 'inches', 'Position', [25 8 6.5 6])

for s = 1:20
    axes('Units', 'inches', 'Position', [pos(s, :) 0.96 1])
    plot(5:10:100, percent_fired(irs_shown(s), :) .* 100, 'k')
    set(gca, 'XTick', [25 50 75])
    
    if ismember(s, [1 6 11 16])
        ylabel('% of cells')
    end
    
    if s >= 16
        xlabel('% replicated')
    end
    
    title(['IR ' num2str(IR_number(irs_shown(s)))])
end

printFigure('out/FigureR6.pdf')
close
clearvars
