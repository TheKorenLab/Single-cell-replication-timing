function replication_state_masked = mask_replication_state(replication_state_filtered, ...
    replication_tracks) 

    load('data/hg37_genome_metadata.mat', 'genome_windows')
    
    replication_state_masked = cell(22, 1);
    for Chr = 1:22
        replication_state_masked{Chr} = replication_state_filtered{Chr};
        barcodes = unique(replication_tracks{Chr}(~isnan(replication_tracks{Chr}(:, 6)), 1));
        
        for barcode = 1:size(replication_state_filtered{Chr}, 2)
            if ~ismember(barcode, barcodes)
                replication_state_masked{Chr}(:, barcode) = NaN;
            end
        end
        
        replication_state_masked{Chr}(replication_state_masked{Chr} == 4) = NaN;
        for track = 1:size(replication_tracks{Chr}, 1)
            
            if isnan(replication_tracks{Chr}(track, 6))
                continue
            end
            
            barcode = replication_tracks{Chr}(track, 1);
            X = [find(genome_windows{Chr}(:, 1) == replication_tracks{Chr}(track, 2)) ...
                find(genome_windows{Chr}(:, 2) == replication_tracks{Chr}(track, 3))];
            replication_state_masked{Chr}(X(1):X(2), barcode) = 4;
        end
    end
 
end
