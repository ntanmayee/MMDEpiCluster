#!/bin/bash

export init_wd=$(pwd)
export SRR_ACC_LIST=SRR_Acc_List.txt
export n_threads=$2

echo $1
cd $1

mkdir fastq
mkdir trimmed
mkdir cleaned
mkdir merged
mkdir tag_dir

less $SRR_ACC_LIST | parallel -j $n_threads "fastq-dump --outdir fastq  --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}"

mv fastq/*.fastq.gz trimmed/

cd trimmed
ls *1.fastq.gz | grep "^[^_]*" -o | parallel -j $n_threads "bowtie2 -x $INDEX -1 {}_pass_1.fastq.gz -2 {}_pass_2.fastq.gz -S {}.sam"
ls *.sam | parallel -j $n_threads "samtools view -bS {} | samtools sort - -o {}.bam"
ls *.bam | parallel -j $n_threads "picard MarkDuplicates \
      I={} \
      O=../cleaned/{} \
      M=../cleaned/marked_dup_metrics_{}.txt \
      REMOVE_DUPLICATES=true"
cd ..

cd merged
samtools merge all.bam $(ls ../cleaned/*.bam)
ls *.bam | parallel -j $n_threads "samtools index {}"
cd ..

makeTagDirectory ./tag_dir $(ls merged/*.bam)
findPeaks tag_dir -style histone -o auto

cd $init_wd