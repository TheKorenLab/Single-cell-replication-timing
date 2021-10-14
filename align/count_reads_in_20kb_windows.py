#! /bin/python

import numpy as np
import pandas as pd
import pysam
import h5py
import argparse

PROFILE_MAPQ_THRESHOLD = 30
PROFILE_INSERT_THRESHOLD = 20000
PROCESSED_BARCODE_TAG = 'CB'
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
    parser.add_argument('--cellranger', dest='cellranger', action='store_true')
    parser.set_defaults(cellranger=False)
    parser.add_argument('--mm10', dest='mm10', action='store_true')
    parser.set_defaults(mm10=False)
    args = parser.parse_args()

    print('Sample: ' + args.id)

    bam_file = pysam.Samfile(args.dir + args.id + '.bam', 'rb')
    assert bam_file.has_index()

    if not args.mm10:
        contigs = [str(chrom) for chrom in range(1, 23)]
        prefix = ''
    else:
        contigs = [str(chrom) for chrom in range(1, 20)]
        prefix = 'chr'

    contigs.extend(['X', 'Y'])

    if not args.cellranger:
        bc_list = pd.read_csv(args.dir + 'barcode_list.csv', header=None)
        bc_list = list(bc_list[0])
    else:
        bc_list = pd.read_csv(args.dir + args.id + '_per_cell_summary_metrics.csv', usecols = [0])
        bc_list = list(bc_list['barcode'])

    nbcs = len(bc_list)
    bc_index = {}
    for bci in range(nbcs):
        bc_index[bc_list[bci]] = bci

    with h5py.File(args.dir + args.id + '.h5', 'w') as store:
        for chrom in contigs: 

            print('Chromosome ' + chrom)
            chrom_length = bam_file.get_reference_length(prefix + chrom)
            nbins = chrom_length//WINDOW_SIZE + int(chrom_length % WINDOW_SIZE !=0)
            counts = np.zeros( (nbcs, nbins), dtype='float32' )
            reads = []

            for rec in bam_file.fetch( prefix + chrom, 0, chrom_length ):

                if not rec.has_tag(PROCESSED_BARCODE_TAG):
                    continue
		
                inserts = inserts_per_alignment(rec)
	
                if inserts == 0:
                    continue
		
                bc = rec.get_tag(PROCESSED_BARCODE_TAG)
                bci = bc_index.get(bc, None)
		
                if bci is None:
                    continue
  
                win = (rec.reference_start+1)//WINDOW_SIZE
                counts[bci, win] += inserts

                reads.append([bci, rec.reference_start+1, inserts])
            
            reads = pd.DataFrame(reads).transpose()
            
            store.create_dataset('reads/' + prefix + chrom, reads.shape, data = reads,
                                 compression='gzip', compression_opts=9)
            store.create_dataset('raw_counts/' + prefix + chrom, counts.shape, data = counts)

if __name__ == '__main__':
    main()
