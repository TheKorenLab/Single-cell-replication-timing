function replication_tracks = identify_replication_tracks(replication_state_filtered)

    replication_tracks = cell(22, 1);
    for Chr = 1:22
        replication_tracks{Chr} = list_cn_segments(replication_state_filtered{Chr}, 4);
        replication_tracks{Chr}(:, 5) = mean(replication_tracks{Chr}(:, 2:3), 2);
        replication_tracks{Chr}(:, 6) = replication_tracks{Chr}(:, 3) - ...
            replication_tracks{Chr}(:, 2) + 1;
    end
    
end
