#!/bin/bash

conda activate fastqProcessing

export fqDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/fastq"
export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bowtie2"
export index="/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly"
export samples="A_H2AZ_H_r1_R1 A_H2AZ_H_r2_R1 A_Inp_H_r1_R1 A_Inp_H_r2_R1 A_shH2AZ_Inp_H_r1_R1 A_shH2AZ_Inp_H_r2_R1 A_TGFb_H2AZ_H_r1_R1 A_TGFb_H2AZ_H_r2_R1 A_TGFb_Inp_H_r1_R1 A_TGFb_Inp_H_r2_R1 CA1a_H2AZ_H_r1_R1 CA1a_H2AZ_H_r2_R1 CA1a_Inp_H_r1_R1 CA1a_Inp_H_r2_R1"
export deeptoolsDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/deeptools"

if [ ! -d $bamDir ]
then
    mkdir -p $bamDir
fi


for sample in $samples; do
    if [ ! -f ${bamDir}/${sample}.bam ]; then
        bowtie2 -x $index -1 ${fqDir}/${sample}.end1.fastq.gz -2 ${fqDir}/${sample}.end2.fastq.gz --threads 28 -X 2000 --no-mixed --no-discordant --no-dovetail --no-contain | samtools view -Sb - > ${bamDir}/${sample}.bam
    fi
done

for sample in $samples; do
    if [ ! -f ${bamDir}/${sample}.sorted.bam ]; then
        samtools sort -@ 32 -m 2G ${bamDir}/${sample}.bam -T ${bamDir}/${sample}.temp -o ${bamDir}/${sample}.sorted.bam
    fi
done

conda deactivate
conda activate bamProcessing

for sample in $samples; do
    if [ ! -f ${bamDir}/${sample}.dedup.bam ]; then
         picard -Djava.io.tmpdir=/tmp MarkDuplicates\
            -INPUT ${bamDir}/${sample}.sorted.bam\
            -OUTPUT ${bamDir}/${sample}.dedup.bam\
            -ASSUME_SORTED TRUE\
            -REMOVE_DUPLICATES TRUE\
            -CREATE_INDEX TRUE\
            -METRICS_FILE ${bamDir}.picard_metrics.txt
    fi
done

for sample in $samples; do
    if [ ! -f ${bamDir}/${sample}.dedup.insert_size_metrics.txt ]; then
        picard CollectInsertSizeMetrics\
            -I ${bamDir}/${sample}.dedup.bam\
            -O ${bamDir}/${sample}.dedup.insert_size_metrics.txt\
            -H ${bamDir}/${sample}.dedup.insert_size_metrics.pdf\
            -AS T
    fi
done

