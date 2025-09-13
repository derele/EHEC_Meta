#!/usr/bin/env bash
# Usage: ./map_sra_bowtie2.sh <SRA_ID>

set -euo pipefail

SRA_ID="$1"
THREADS=4
REF_PREFIX="VFDB_setA_nt_stx_eae"   # bowtie2 reference prefix (bt2 index files)
DB="VFDB_setA_pro_stx_eaeDB.dmnd"   # DIAMOND protein DB

### DOWNLOAD
echo "▶ Downloading $SRA_ID via prefetch..."
prefetch "$SRA_ID"

echo "▶ Converting to FASTQ..."
fasterq-dump --threads 8 --skip-technical --split-files ${SRA_ID}/${SRA_ID}.sra

### CHECK IF SINGLE-END OR PAIRED-END
if [[ -f ${SRA_ID}_1.fastq && -f ${SRA_ID}_2.fastq ]]; then
    echo "▶ Detected paired-end data"

    repair.sh in1=${SRA_ID}_1.fastq in2=${SRA_ID}_2.fastq \
              out1=${SRA_ID}_1.repaired.fastq \
              out2=${SRA_ID}_2.repaired.fastq \
              outs=${SRA_ID}_singletons.fastq \
              overwrite=t

    echo "▶ Running fastp..."
    fastp -i ${SRA_ID}_1.repaired.fastq -I ${SRA_ID}_2.repaired.fastq \
          -o ${SRA_ID}_1.clean.fastq -O ${SRA_ID}_2.clean.fastq \
          --thread $THREADS

    READ1=${SRA_ID}_1.clean.fastq
    READ2=${SRA_ID}_2.clean.fastq
else
    echo "▶ Detected single-end data"

    # manchmal heißt die Datei ${SRA_ID}.fastq
    if [[ -f ${SRA_ID}.fastq ]]; then
        mv ${SRA_ID}.fastq ${SRA_ID}_SE.fastq
    fi

    echo "▶ Running fastp..."
    fastp -i ${SRA_ID}_SE.fastq -o ${SRA_ID}_SE.clean.fastq \
          --thread $THREADS

    READ1=${SRA_ID}_SE.clean.fastq
    READ2=""
fi

### BOWTIE2
echo "▶ Mapping with bowtie2..."
if [[ -n "$READ2" ]]; then
    bowtie2 -x "$REF_PREFIX" \
            -1 $READ1 -2 $READ2 \
            --very-sensitive -p $THREADS | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o ${SRA_ID}.bt2.bam
else
    bowtie2 -x "$REF_PREFIX" \
            -U $READ1 \
            --very-sensitive -p $THREADS | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o ${SRA_ID}.bt2.bam
fi

echo "▶ Indexing and summarizing..."
samtools index ${SRA_ID}.bt2.bam
samtools idxstats ${SRA_ID}.bt2.bam > ${SRA_ID}.bt2.idxstats.txt

### DIAMOND
echo "▶ Preparing input for DIAMOND..."
if [[ -n "$READ2" ]]; then
    cat $READ1 $READ2 > ${SRA_ID}.clean.fastq
else
    cp $READ1 ${SRA_ID}.clean.fastq
fi

echo "▶ Running DIAMOND blastx..."
diamond blastx -d "$DB" \
        -q ${SRA_ID}.clean.fastq \
        --out ${SRA_ID}.diamond.tsv \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --evalue 1e-10 --max-target-seqs 1 --threads $THREADS

echo "✅ Finished processing $SRA_ID"

### CLEANUP
rm -f ${SRA_ID}_*.fastq
rm -f ${SRA_ID}.bt2.bam*
rm -rf ${SRA_ID}
rm -f ${SRA_ID}.clean.fastq
rm -f fastp.html fastp.json
