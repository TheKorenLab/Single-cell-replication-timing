#! /usr/bin/env python3

import sys
import gzip
import argparse
import numpy as np
from ctypes import c_int
from multiprocessing import Pool, Value, Lock, Semaphore
from Bio.SeqIO.QualityIO import FastqGeneralIterator

VALID_BARCODES = None

def init(valid_barcodes):
    global VALID_BARCODES
    VALID_BARCODES = valid_barcodes

def convert_phred_to_qscore(quality_score):
    return np.array([ord(c)-33 for c in quality_score])

def list_hamming_cousins(raw_barcode, raw_quality_score):
    index = np.where(raw_quality_score < 24)[0][0]
    return {raw_barcode[:index] + nt + raw_barcode[index+1:] for nt in 'ACGT'}

def select_best_barcode(raw_barcode, raw_quality_score):
    if raw_barcode in VALID_BARCODES:
        return raw_barcode
    if sum(raw_quality_score < 24) != 1:
        return None
    cousins = list_hamming_cousins(raw_barcode, raw_quality_score)
    possible_barcodes = cousins.intersection(VALID_BARCODES)
    if len(possible_barcodes) == 1:
        return possible_barcodes.pop()
    return None

def reader(fastq_R1, fastq_R2, generator_buffer):
   for (name_R1, seq_R1, qual_R1), (name_R2, seq_R2, qual_R2) \
         in zip(FastqGeneralIterator(fastq_R1),FastqGeneralIterator(fastq_R2)):
         generator_buffer.acquire()
         yield (name_R1, seq_R1, qual_R1, name_R2, seq_R2, qual_R2)    

def process_record(inputs):
    name_R1, seq_R1, qual_R1, name_R2, seq_R2, qual_R2 = inputs
    name_R1 = str.split(name_R1)[0]
    name_R2 = str.split(name_R2)[0]
    assert name_R1 == name_R2
    barcode = select_best_barcode(seq_R1[:16], convert_phred_to_qscore(qual_R1[:16]))
    if barcode is not None:
        return('@%s:%s\n%s\n+\n%s\n@%s:%s\n%s\n+\n%s\n'
                 % (name_R1, barcode, seq_R1[16:], qual_R1[16:], name_R2, barcode, seq_R2, qual_R2))

def return_processed_record(generator_buffer, result):
    generator_buffer.release()
    return result
    
class Counter:
    def __init__(self):
        self.data = Value(c_int)
        self.lock = Lock()     
    
    def increment(self):
        with self.lock:
            self.data.value += 1
    
    def evaluate(self):
        return self.data.value

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--id')
    parser.add_argument('--processes', type=int)
    parser.add_argument('--buffer', type=int)
    args = parser.parse_args()

    num_read_pairs = Counter()
    num_filtered_pairs = Counter()
    
    valid_barcodes = frozenset([x.strip() for x in open('barcode_list.csv', 'r') if not ('#' in x)])

    generator_buffer =  Semaphore(args.buffer) 
    with gzip.open(args.id + '_R1.fastq.gz', 'rt') as read1, gzip.open(args.id + '_R2.fastq.gz', 'rt') as read2:
        with Pool(processes=args.processes, initializer=init, initargs=(valid_barcodes,)) as p:
            for result in p.imap(process_record, reader(read1, read2, generator_buffer), chunksize=10):
                processed_record = return_processed_record(generator_buffer, result)
                num_read_pairs.increment()
                if result is not None:            
                    num_filtered_pairs.increment()
                    print(processed_record, end='')
    
    print('Total read pairs: ' + f"{num_read_pairs.evaluate():,}", file=sys.stderr)
    print('Total filtered read pairs: ' + f"{num_filtered_pairs.evaluate():,}", file=sys.stderr)
    
if __name__ == '__main__':
    main()

