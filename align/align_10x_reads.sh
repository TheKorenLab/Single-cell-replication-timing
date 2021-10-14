#! /bin/bash

set -eo pipefail
source logger.sh

REF_GENOME=$1
ID=$2

log 'Sample: ' $ID
 
# MERGE INPUT FILES
REPS=$(ls $INPUT_PATH"/"$ID*"R1.fastq.gz" | wc -l)

if [ $REPS -gt 1 ]; then
    log "Merging" $REPS "sequencing replicates"
    cat $INPUT_PATH"/"$ID*"R1.fastq.gz" > $ID"_R1.fastq.gz" &
    cat $INPUT_PATH"/"$ID*"R2.fastq.gz" > $ID"_R2.fastq.gz" &
    wait
else
    cp $INPUT_PATH"/"$ID"-replicate1_R1.fastq.gz" $ID'_R1.fastq.gz' &
    cp $INPUT_PATH"/"$ID"-replicate1_R2.fastq.gz" $ID'_R2.fastq.gz' &
    wait
fi

# NAIVE COUNT OF BARCODES
log "Counting barcodes"
$SEQTK trimfq -L16 $ID'_R1.fastq.gz' | \
    $SCRIPTS"/naive_count_barcodes.py" /dev/stdin
log $(wc -l < barcode_list.csv) "cells detected."

# FILTER AND ALIGN READS
log 'Filtering and aligning reads'
$SCRIPTS"/handle_barcodes.py" --id $ID --processes 16 --buffer 2000000 | \
  bwa mem -M -t8 -p $REF_GENOME /dev/stdin > \
  $ID".unsorted.bam" 2> bwa.log

log 'bwa complete'

# SORT BAM FILE
log 'Sorting BAM file'
samtools sort -@8 $ID".unsorted.bam" -o $ID'.bam' 2>> bwa.log
rm $ID".unsorted.bam"

# MERGE BARCODES INTO READS
log 'Merging barcodes'
samtools view -h $ID'.bam' | \
  awk -f $SCRIPTS"/tag_reads_with_barcodes.awk" | \
  samtools view -Shb - > $ID'.barcoded.bam'
rm $ID".bam"
 
# BARCODE-AWARE DUPLICATE MARKING
log 'Marking duplicates per cell'
java -jar $PICARD UmiAwareMarkDuplicatesWithMateCigar \
  INPUT=$ID'.barcoded.bam' \
  OUTPUT=$ID'.barcoded.dups.marked.bam' \
  METRICS_FILE=/dev/null \
  UMI_TAG_NAME='CB' \
  TMP_DIR='/SSD/temp_bwa/' 2> picard.log
mv $ID".barcoded.dups.marked.bam" $ID".bam"
rm $ID".barcoded.bam"

# INDEX BAM FILE
samtools index $ID".bam" $ID".bai"

# CLEAN UP
rm $ID"_R1.fastq.gz" $ID"_R2.fastq.gz"
chmod 777 *
log $0 "complete."

exit 0
