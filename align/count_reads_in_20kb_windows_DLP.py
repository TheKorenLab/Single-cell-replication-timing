#! /bin/python

import numpy as np
import pandas as pd
import pysam
import h5py
import argparse

PROFILE_MAPQ_THRESHOLD = 30
PROFILE_INSERT_THRESHOLD = 20000
WINDOW_SIZE = int(20000)

def inserts_per_alignment(rec):
    if rec.is_secondary or rec.is_supplementary:
        return 0.0
    if (rec.is_unmapped or
        rec.mapping_quality < PROFILE_MAPQ_THRESHOLD or
        rec.is_duplicate):
        return 0.0
    if rec.mate_is_unmapped:
        return 1.0
    if (rec.reference_id == rec.next_reference_id):
        if rec.template_length > PROFILE_INSERT_THRESHOLD:
            return 1.0
        else:
            return 0.5
    else:
        return 1.0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir')
    parser.add_argument('--id')
    args = parser.parse_args()

    bc_list = pd.read_csv(args.dir + 'bam.list', header=None)
    bc_list = list(bc_list[0])
    nbcs = len(bc_list)

    contigs = [str(chrom) for chrom in range(1, 23)]
    contigs.extend(['X', 'Y'])
    prefix = 'chr'

    with h5py.File(args.dir + args.id + '.h5', 'w') as store:
        for chrom in contigs: 
            
            print('Chromosome ' + chrom)
      
            bam_file = pysam.Samfile(args.dir + bc_list[0], 'rb')
            assert bam_file.has_index()
            chrom_length = bam_file.get_reference_length(chrom)
            nbins = chrom_length//WINDOW_SIZE + int(chrom_length % WINDOW_SIZE !=0)
            counts = np.zeros( (nbcs, nbins), dtype='float32' )
            reads = []

            for barcode in range(nbcs):
                bam_file = pysam.Samfile(args.dir + bc_list[barcode], 'rb')
                assert bam_file.has_index()

                for rec in bam_file.fetch( chrom, 0, chrom_length ):

                    inserts = inserts_per_alignment(rec)
	
                    if inserts == 0:
                        continue
		
                    win = (rec.reference_start+1)//WINDOW_SIZE
                    counts[barcode, win] += inserts

                    reads.append([barcode, rec.reference_start+1, inserts])
            
            reads = pd.DataFrame(reads).transpose()
            
            store.create_dataset('reads/' + prefix + chrom, reads.shape, data = reads,
                                 compression='gzip', compression_opts=9)
            store.create_dataset('raw_counts/' + prefix + chrom, counts.shape, data = counts)

if __name__ == '__main__':
    main()
