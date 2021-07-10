#!/bin/sh
# Script to sort TIER-Seq non-deduplicated and unsorted reads by TSS / PSS tags, unassigned: transcript reads

## set working directory
# directory in which data should be saved
WORK_DIR="/data/documents/TIER_Seq/pre_processing/sorted_not-deduplicated/" 
# directory where to find bbduk
BBDUK="/home/utee/bbmap/bbduk.sh" 
# directory in which unsorted data is located
INPUT_DIR="/data/documents/TIER_Seq/pre_processing/not-deduplicated/"

# PSS and TSS tags
PSS="PSS.fa"
TSS="TSS.fa"

# actual sorting
cd $INPUT_DIR
for File in `find -type f -name \*.fastq.gz`; do 
# run bbduk to filter out PSS
bash $BBDUK in=$File out=$WORK_DIR"unassigned/"${File%.fastq.gz}"_unassigned_noPSS.fastq.gz" outm=$WORK_DIR"PSS/"${File%.fastq.gz}"_PSS_noTrim.fastq.gz" stats=$WORK_DIR"stats/"${File%.fastq.gz}"_PSS_stats.txt" ref=$WORK_DIR$PSS statscolumns=5 restrictleft=23 k=8

# run bbduk to filter out TSS
bash $BBDUK in=$WORK_DIR"unassigned/"${File%.fastq.gz}"_unassigned_noPSS.fastq.gz" out=$WORK_DIR"unassigned/"${File%.fastq.gz}"_unassigned.fastq.gz" outm=$WORK_DIR"TSS/"${File%.fastq.gz}"_TSS_noTrim.fastq.gz" stats=$WORK_DIR"stats/"${File%.fastq.gz}"_TSS_stats.txt" ref=$WORK_DIR$TSS statscolumns=5 restrictleft=23 k=8

rm $WORK_DIR"unassigned/"${File%.fastq.gz}"_unassigned_noPSS.fastq.gz"

# filter out PSS tags and 2 nt
bash $BBDUK in=$WORK_DIR"PSS/"${File%.fastq.gz}"_PSS_noTrim.fastq.gz" out=$WORK_DIR"PSS/"${File%.fastq.gz}"_PSS.fastq.gz" stats=$WORK_DIR"stats/"${File%.fastq.gz}"_PSS_trimming_stats.txt" ref=$WORK_DIR$PSS statscolumns=5 restrictleft=23 k=8 ktrim=l tp=2 overwrite=t

rm $WORK_DIR"PSS/"${File%.fastq.gz}"_PSS_noTrim.fastq.gz"

# filter out TSS tags and 2 nt
bash $BBDUK in=$WORK_DIR"TSS/"${File%.fastq.gz}"_TSS_noTrim.fastq.gz" out=$WORK_DIR"TSS/"${File%.fastq.gz}"_TSS.fastq.gz" stats=$WORK_DIR"stats/"${File%.fastq.gz}"_TSS_trimming_stats.txt" ref=$WORK_DIR$TSS statscolumns=5 restrictleft=23 k=8 ktrim=l tp=2 overwrite=t

rm $WORK_DIR"TSS/"${File%.fastq.gz}"_TSS_noTrim.fastq.gz"
done



