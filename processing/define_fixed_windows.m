% Define 20kb fixed-size windows for hg37

% Chromosome sizes downloaded from
% http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
fid = fopen('data/hg37.chrom.sizes.txt', 'r');
infile = textscan(fid, '%s%f', 'Delimiter', '\t');
fclose(fid);

chrom_size = NaN(22, 1);
for Chr = 1:22
    is_chr = strcmp(infile{1}, ['chr' num2str(Chr)]);
    chrom_size(Chr) = infile{2}(is_chr);
end

genome_windows = cell(22, 1);
for Chr = 1:22
    genome_windows{Chr} = NaN(ceil(chrom_size(Chr) / 20e3), 3);
    genome_windows{Chr}(1:end-1, 1) = 0:20e3:chrom_size(Chr)-20e3;
    genome_windows{Chr}(1:end-1, 2) = 20e3:20e3:chrom_size(Chr);
    genome_windows{Chr}(end, 1:2) = [genome_windows{Chr}(end-1, 2) chrom_size(Chr)];
    genome_windows{Chr}(:, 3) = round(mean(genome_windows{Chr}(:, 1:2), 2));
end

save('data/hg37_genome_metadata.mat', 'chrom_size', 'genome_windows')
clearvars -except samples
