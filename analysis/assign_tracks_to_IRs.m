function [replication_tracks, single_cell_IRs] = assign_tracks_to_IRs(replication_tracks, ...
    single_cell_IRs)

    for Chr = 1:22

        cluster_id = NaN(size(replication_tracks{Chr}, 1), 1);

        for id = 1:size(single_cell_IRs{Chr}, 1)
           is_in_cluster = replication_tracks{Chr}(:, 4) >= ...
               single_cell_IRs{Chr}(id, 2) & ...
               replication_tracks{Chr}(:, 4) <= ...
               single_cell_IRs{Chr}(id, 3) & ...
               replication_tracks{Chr}(:, 5) < 1e6;
           cluster_id(is_in_cluster) = id;
       end
       replication_tracks{Chr}(:, 8) = cluster_id;

       defines_IR = ~isnan(replication_tracks{Chr}(:, 6));
       is_mismatched = replication_tracks{Chr}(defines_IR, 6) ~= ...
           replication_tracks{Chr}(defines_IR, 8) & ...
           ~isnan(replication_tracks{Chr}(defines_IR, 8));
       assert(~any(is_mismatched))

       replication_tracks{Chr}(defines_IR, 8) = replication_tracks{Chr}(defines_IR, 6);

       single_cell_IRs{Chr}(:, 12) = NaN;
       for id = 1:size(single_cell_IRs{Chr}, 1)
           single_cell_IRs{Chr}(id, 12) = sum(replication_tracks{Chr}(:, 8) == id);
       end

       assert(all(single_cell_IRs{Chr}(:, 12) >= single_cell_IRs{Chr}(:, 10)))

    end
    
end

