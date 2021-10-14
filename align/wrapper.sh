#! /bin/bash

usage() {
echo "Usage: --ref PATH/TO/REFERENCE/GENOME"
    echo "       --id  PATH/TO/INPUT/FILE"
    exit 1
}

# PARSE ARGUMENTS
while [ "$1" != "" ]; do
    case $1 in
        --ref )    shift
                   REF_GENOME=$1
                   ;;
       --id  )     shift
                   INPUT_FILE=$1
                   ;;
         *   )     usage
                   ;;
    esac
    shift
done

# VALIDATE INPUTS
if [ -z $REF_GENOME ] || [ -z $INPUT_FILE ]; then
    usage
fi

if [ ! -f $REF_GENOME'.fa' ] || [ ! -f $INPUT_FILE ]; then
    echo "One or more input files is not valid."
    exit 1
fi

# EXTRACT SAMPLE NAME
SAMPLE="${INPUT_FILE##*/}"
ID="${SAMPLE%-replicate1_R1*}"

# SET OUTDIR
OUTDIR="./"$ID
if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi
cd $OUTDIR

# ALIGN READS TO REFERENCE GENOME
export SCRIPTS="align"
export SEQTK="/PATH/TO/SEQTK"
export PICARD="/PATH/TO/PICARDTOOLS"
export INPUT_PATH=$(dirname $INPUT_FILE)

export LOGFILE=log
> $LOGFILE

./align_10x_reads.sh $REF_GENOME $ID > log 2>&1
STATUS=$?
