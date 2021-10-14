function [p_cells_fired, late_IR_class] = classify_late_IRs(single_cell_IRs, ...
    replication_state_filtered, percent_replicated_filtered, replication_tracks)

    load('data/hg37_genome_metadata.mat', 'genome_windows')
    
    p_cells_fired = cell(22, 1);
    late_IR_class = cell(22, 1);
    for Chr = 1:22
        
        % Pull out windows containing IR centers
        IR_centers = NaN(size(single_cell_IRs{Chr}, 1), 1);
        for o = 1:size(single_cell_IRs{Chr}, 1)
            IR_centers(o) = find(genome_windows{Chr}(:, 1) < single_cell_IRs{Chr}(o, 4) & ...
                genome_windows{Chr}(:, 2) >= single_cell_IRs{Chr}(o, 4));
        end
        
        late_IR_class{Chr} = false(size(IR_centers, 1), 3);
        p_cells_fired{Chr} = NaN(size(IR_centers, 1), 1);
        
        for o = 1:size(IR_centers, 1)

            if single_cell_IRs{Chr}(o, 10) < 2
                continue
            end

            IR_state = replication_state_filtered{Chr}(IR_centers(o), :);
            p_cells_fired{Chr}(o) = sum(IR_state == 4) / sum(~isnan(IR_state));
            
            if p_cells_fired{Chr}(o) > 0.5
                continue
            end        

            rep_cells = replication_tracks{Chr}(replication_tracks{Chr}(:, 8) == o, 1);
            p_rep_cells = percent_replicated_filtered(rep_cells);

            late_IR_class{Chr}(o, 1) = all(p_rep_cells > 0.5);
            late_IR_class{Chr}(o, 2) = sum(p_rep_cells < 0.5) > 0 & sum(p_rep_cells < 0.5) <= 5;
            late_IR_class{Chr}(o, 3) = sum(p_rep_cells < 0.5) > 5 & ~all(p_rep_cells > 0.5);
        end        
    end

end
