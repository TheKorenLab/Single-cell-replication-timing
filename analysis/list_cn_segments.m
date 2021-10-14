function replicated_segments = list_cn_segments(replication_state_matrix, cn)

    num_cells = size(replication_state_matrix, 2);
    
    is_gap = sum(isnan(replication_state_matrix), 2) == num_cells;
    
    replicated_segments = cell(num_cells, 1);
    for barcode = 1:num_cells

        values = replication_state_matrix(:, barcode);
        values(:, 2) = 1:length(values);

        values(:, 3) = 0;
        values(1, 3) = 1;

        v = values(1,1);
        na_counter = 0;
        for row = 1:size(values, 1)

            if is_gap(row)
                values(row, 3) = 1;
                continue
            end
            
            if isnan(values(row, 1))
                na_counter = na_counter + 1;

                if na_counter > 5
                    values(row, 3) = 1;
                    na_counter = 0;
                end
                continue
            else
                na_counter = 0;
            end

            if values(row, 1) ~= v
                values(row, 3) = 1;
                v = values(row, 1);
            end
        end
                
        segments = values(values(:, 3) == 1, 2);
        segments(1:end-1, 2) = segments(2:end, 1) - 1;
        segments(end, 2) = values(end, 2);
               
        segments(:, 3) = NaN;
        for row = 1:size(segments, 1)
            segment_cn = values(segments(row, 1):segments(row, 2), 1);
            segment_cn = unique(segment_cn(~isnan(segment_cn)));

            if ~isempty(segment_cn)
                segments(row, 3) = segment_cn;
            end
        end

        % Flag rows that are neighbored by an NaN
        segments(:, 4) = 0;
        for row = 1:size(segments, 1)
            neighbors = [row-1 row+1];
            neighbors = neighbors(neighbors > 0 & neighbors < size(segments, 1)-1);

            segments(row, 4) = ~isnan(segments(row, 3)) & ...
                any(isnan(segments(neighbors, 3)));
        end

        for row = 1:size(segments,1)-1
            if segments(row, 3) == 4 && segments(row+1, 3) == 4
                segments(row:row+1, 4) = 1;
            end
        end
        
        % Flag rows that contain too many NaNs
        segments(:, 5) = segments(:, 2) - segments(:, 1) + 1;
        segments(:, 6) = 0;
        for row = 1:size(segments, 1)
            segments(row, 6) = sum(isnan(values(segments(row, 1):segments(row, 2), 1)));
        end
        segments(:, 6) = segments(:, 6) ./ segments(:, 5);
        has_too_many_NaNs = segments(:, 6) >= 0.5;
        segments(has_too_many_NaNs, 4) = 1;
        
        % Return only the segments of the desired copy number
        is_cn4 = segments(:, 3) == cn;
        
        replicated_segments{barcode}(:, 2:4) = segments(is_cn4, [1 2 4]);
        replicated_segments{barcode}(:, 1) = barcode;
       
    end

    replicated_segments = cell2mat(replicated_segments);

end
