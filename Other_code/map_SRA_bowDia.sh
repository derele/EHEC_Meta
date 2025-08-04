#!/bin/bash
# Usage: ./map_sra_bowtie2.sh <SRA_ID>

set -euo pipefail

SRA_ID="$1"
THREADS=4
REF_PREFIX="VFDB_setA_nt_stx_eae" # bowtie2 reference prefix (bt2 index files)
DB="VFDB_setA_pro_stx_eaeDB.dmnd" # DIAMOND protein DB

### DOWNLOAD

echo "▶ Downloading $SRA_ID via prefetch..."
prefetch "$SRA_ID"

echo "▶ reading from ${SRA_ID}/${SRA_ID}.sra"
# 2. Convert in RAM or SSD
fasterq-dump --threads 8 --skip-technical --split-files \
  ${SRA_ID}/${SRA_ID}.sra

### CLEANUP

repair.sh in1=${SRA_ID}_1.fastq in2=${SRA_ID}_2.fastq \
          out1=${SRA_ID}_1.repaired.fastq \
          out2=${SRA_ID}_2.repaired.fastq \
          outs=${SRA_ID}_singletons.fastq \
                    overwrite=t


echo "▶ Running fastp..."
fastp -i ${SRA_ID}_1.repaired.fastq -I ${SRA_ID}_2.repaired.fastq \
      -o ${SRA_ID}_1.clean.fastq -O ${SRA_ID}_2.clean.fastq \
      --thread $THREADS 

### BOWTIE

echo "▶ Mapping with bowtie2..."
bowtie2 -x "$REF_PREFIX" \
        -1 ${SRA_ID}_1.clean.fastq -2 ${SRA_ID}_2.clean.fastq \
        --very-sensitive -p $THREADS | \
    samtools view -bS - | \
    samtools sort -@ $THREADS -o ${SRA_ID}.bt2.bam

echo "▶ Indexing and summarizing..."
samtools index ${SRA_ID}.bt2.bam
samtools idxstats ${SRA_ID}.bt2.bam > ${SRA_ID}.bt2.idxstats.txt

### DIAMOND
echo "▶ Preparing single FASTQ for DIAMOND..."
cat ${SRA_ID}_1.clean.fastq ${SRA_ID}_2.clean.fastq > ${SRA_ID}.clean.fastq

echo "▶ Running DIAMOND blastx..."
diamond blastx -d "$DB" \
        -q ${SRA_ID}.clean.fastq \
        --out ${SRA_ID}.diamond.tsv \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --evalue 1e-10 --max-target-seqs 1 --threads $THREADS

echo "✅ Finished processing $SRA_ID"

# Optional cleanup
rm ${SRA_ID}_*.fastq
rm ${SRA_ID}.bt2.bam*

rm -r ${SRA_ID}
rm ${SRA_ID}.clean.fastq

rm fastp.html fastp.json

