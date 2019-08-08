#!/bin/bash

mkdir fastq
mkdir report
mkdir trimmed
mkdir annotation
mkdir cleaned
mkdir tag_dir

export SRR_ACC_LIST=SRR_Acc_List.txt

less $SRR_ACC_LIST | parallel -j 45 "fastq-dump --outdir fastq  --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}"

ls fastq/*fastq.gz | parallel -j 45 "trim_galore  --stringency 3 --fastqc -o trimmed/ {}"

mv trimmed/*.html report/
mv trimmed/*.txt report/
mv trimmed/*.zip report/

multiqc report

cd annotation
wget ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna//Mus_musculus.GRCm38.dna.toplevel.fa.gz
cd ..

cd trimmed
bowtie2-build ../annotation/Mus_musculus.GRCm38.dna.toplevel.fa.gz mm9
bowtie2-inspect -n mm9
cd ..

cd trimmed
ls *_trimmed.fq.gz | parallel -j 45 "bowtie2 -x mm9 -U {} -S {}.sam"
ls *.sam | parallel -j 45 "samtools view -bS {} | samtools sort - -o {}.bam"
ls *.bam | parallel -j 45 "picard MarkDuplicates \
      I={} \
      O=../cleaned/{}.bam \
      M=../cleaned/marked_dup_metrics_{}.txt \
      REMOVE_DUPLICATES=true"
cd ..

cd cleaned
ls *.bam | parallel -j 45 "samtools index {}"
cd ..

ls cleaned/*bam | makeTagDirectory ./tag_dir/ {}
findPeaks Peaks-ChIP-Seq/ -style histone -o auto -i Input-ChIP-Seq