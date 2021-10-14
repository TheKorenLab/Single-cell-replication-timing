function count_reads_in_G1_windows(hdf_file, sample_name)

    load('data/hg37_genome_metadata.mat', 'is_mappable')

    coverage_20kb = cell(22, 1);
    for Chr = 1:22
        coverage_20kb{Chr} = double(h5read(hdf_file, ['/raw_counts/chr' num2str(Chr)]));
        coverage_20kb{Chr}(~is_mappable{Chr}, :) = NaN;
    end
    
    [mean_coverage_1Mb, scaled_mapd_1Mb, is_G1, g1_barcodes] = identify_G1_cells(coverage_20kb, ...
        is_mappable);
    num_cells = length(mean_coverage_1Mb);
    
   [g1_windows, is_masked_window, raw_counts_20kb] = define_G1_windows(hdf_file, g1_barcodes, ...
       num_cells);
   
   save('-v7.3', ['data/intermediate/' sample_name '.mat'], 'mean_coverage_1Mb', 'scaled_mapd_1Mb', ...
       'is_G1', 'g1_barcodes', 'g1_windows', 'is_masked_window', 'raw_counts_20kb')

end
