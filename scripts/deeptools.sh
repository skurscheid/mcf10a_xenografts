#!/bin/bash

conda activate deeptools

export fqDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/fastq"
export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bowtie2"
export index="/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly"
export samples="A_H2AZ_H_r1_R1 A_H2AZ_H_r2_R1 A_Inp_H_r1_R1 A_Inp_H_r2_R1 A_shH2AZ_Inp_H_r1_R1 A_shH2AZ_Inp_H_r2_R1 A_TGFb_H2AZ_H_r1_R1 A_TGFb_H2AZ_H_r2_R1 A_TGFb_Inp_H_r1_R1 A_TGFb_Inp_H_r2_R1 CA1a_H2AZ_H_r1_R1 CA1a_H2AZ_H_r2_R1 CA1a_Inp_H_r1_R1 CA1a_Inp_H_r2_R1"
export deeptoolsDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/deeptools"

if [ ! -d $deeptoolsDir/bamCoverage ]; then
    mkdir -p $deeptoolsDir/bamCoverage
fi

for sample in $samples; do
    if [ ! -f ${deeptoolsDir}/bamCoverage/${sample}.bw ]; then
        bamCoverage \
            --effectiveGenomeSize 40999507\
            --binSize 1\
            --centerReads\
            --numberOfProcessors 32\
            --MNase\
            --bam ${bamDir}/${sample}.dedup.bam\
            --outFileName ${deeptoolsDir}/bamCoverage/${sample}.bw
    fi
done

if [ ! -d ${deeptoolsDir}/bigwigCompare ]; then
    mkdir -p ${deeptoolsDir}/bigwigCompare
fi

bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bamCoverage/CA1a_H2AZ_H_r1_R1.bw\
    --bigwig2 ${deeptoolsDir}/bamCoverage/CA1a_Inp_H_r1_R1.bw\
    --skipZeroOverZero\
    --skipNonCoveredRegions\
    --binSize 10\
    --numberOfProcessors 32\
    --outFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r1_R1.bw

bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bamCoverage/CA1a_H2AZ_H_r2_R1.bw\
    --bigwig2 ${deeptoolsDir}/bamCoverage/CA1a_Inp_H_r2_R1.bw\
    --skipZeroOverZero\
    --skipNonCoveredRegions\
    --binSize 10\
    --numberOfProcessors 32\
    --outFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r2_R1.bw

  