#!/bin/bash
export fqDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/fastq"
export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bowtie2"
export index="/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly"
export deeptoolsDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/deeptools"
export regionsDir="/home/sebastian/Data/Tremethick/Breast/Xenografts"
export rootDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap"

#export samples="A_H2AZ_H_r1_R1 A_H2AZ_H_r2_R1 A_Inp_H_r1_R1 A_Inp_H_r2_R1 A_shH2AZ_Inp_H_r1_R1 A_shH2AZ_Inp_H_r2_R1 A_TGFb_H2AZ_H_r1_R1 A_TGFb_H2AZ_H_r2_R1 A_TGFb_Inp_H_r1_R1 A_TGFb_Inp_H_r2_R1 CA1a_H2AZ_H_r1_R1 CA1a_H2AZ_H_r2_R1 CA1a_Inp_H_r1_R1 CA1a_Inp_H_r2_R1"
export samples="CA1a_H2AZ_H_r1_R1 CA1a_H2AZ_H_r2_R1 CA1a_Inp_H_r1_R1 CA1a_Inp_H_r2_R1"

conda activate deeptools

if [ ! -d $deeptoolsDir/bamCoverage ]; then
    mkdir -p $deeptoolsDir/bamCoverage
fi

export samples="A_shH2AZ_Inp_H_r1_R1 A_shH2AZ_Inp_H_r2_R1"
for sample in $samples; do
   # if [ ! -f ${deeptoolsDir}/bamCoverage/${sample}.bw ]; then
        bamCoverage \
            --effectiveGenomeSize 40999507\
            --binSize 10\
            --smoothLength 15\
            --centerReads\
            --numberOfProcessors 32\
            --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
            --bam ${bamDir}/${sample}.dedup.bam\
            --outFileName ${deeptoolsDir}/bamCoverage/${sample}_noncovered.bw
   # fi
done

if [ ! -d ${deeptoolsDir}/bigwigCompare ]; then
    mkdir -p ${deeptoolsDir}/bigwigCompare
fi
# consolidate replicates into merged bigwigs for plotting
export samples="A_H2AZ_H A_Inp_H A_shH2AZ_Inp_H A_TGFb_H2AZ_H A_TGFb_Inp_H CA1a_H2AZ_H CA1a_Inp_H"
export samples="A_shH2AZ_Inp_H"

for i in $samples; do ls ${deeptoolsDir}/bamCoverage/${sample}*_noncovered.bw; done

for sample in $samples; do
    if [ -f ${deeptoolsDir}/bamCoverage/${sample}_r1_R1_noncovered.bw ];
    then
        bigwigCompare \
            --bigwig1 ${deeptoolsDir}/bamCoverage/${sample}_r1_R1_noncovered.bw\
            --bigwig2 ${deeptoolsDir}/bamCoverage/${sample}_r2_R1_noncovered.bw\
            --skipZeroOverZero \
            --operation mean\
            --binSize 10\
            --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
            --numberOfProcessors 32\
            --outFileName ${deeptoolsDir}/bigwigCompare/${sample}_noncovered.bw
    else 
        echo 'input bigwig for ' ${sample} ' missing'
    fi
done

# run bigwigCompare on pooled samples
bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_noncovered.bw\
    --bigwig2 ${deeptoolsDir}/bigwigCompare/CA1a_Inp_H_noncovered.bw\
    --binSize 10\
    --numberOfProcessors 32\
    --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
    --outFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_vs_Inp_H.bw

bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bigwigCompare/A_TGFb_H2AZ_H_noncovered.bw\
    --bigwig2 ${deeptoolsDir}/bigwigCompare/A_TGFb_Inp_H_noncovered.bw\
    --binSize 10\
    --numberOfProcessors 32\
    --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
    --outFileName ${deeptoolsDir}/bigwigCompare/A_TGFb_H2AZ_vs_Inp_H.bw

bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bigwigCompare/A_H2AZ_H_noncovered.bw\
    --bigwig2 ${deeptoolsDir}/bigwigCompare/A_Inp_H_noncovered.bw\
    --binSize 10\
    --numberOfProcessors 32\
    --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
    --outFileName ${deeptoolsDir}/bigwigCompare/A_H2AZ_vs_Inp_H.bw
# 
bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bamCoverage/CA1a_H2AZ_H_r1_R1_noncovered.bw\
    --bigwig2 ${deeptoolsDir}/bamCoverage/CA1a_Inp_H_r1_R1_noncovered.bw\
    --binSize 10\
    --numberOfProcessors 32\
    --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
    --outFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r1_R1.bw

bigwigCompare \
    --bigwig1 ${deeptoolsDir}/bamCoverage/CA1a_H2AZ_H_r2_R1_noncovered.bw\
    --bigwig2 ${deeptoolsDir}/bamCoverage/CA1a_Inp_H_r2_R1_noncovered.bw\
    --binSize 10\
    --numberOfProcessors 32\
    --blackListFileName ${rootDir}/non_covered_hg38_ensembl.bed\
    --outFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r2_R1.bw

export regions="LZ_vs_WT SZ_vs_WT"

for region in $regions; do

    if [ ! -d ${deeptoolsDir}/plotHeatmap/TSS/${region}/ ]; then
        mkdir -p ${deeptoolsDir}/plotHeatmap/TSS/${region}/
    fi

    if [ ! -d ${deeptoolsDir}/computeMatrix/TSS/${region}/ ]; then
        mkdir -p ${deeptoolsDir}/computeMatrix/TSS/${region}/
    fi

    computeMatrix reference-point \
        --regionsFileName ${regionsDir}/${region}_cluster1.bed ${regionsDir}/${region}_cluster2.bed\
        --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r1_R1.bw\
        --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
        --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.bed\
        --referencePoint TSS\
        --beforeRegionStartLength 1000\
        --afterRegionStartLength 1000\
        --binSize 10\
        --missingDataAsZero\
        --sortRegions keep

    plotHeatmap \
        --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
        --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r1_R1.pdf\
        --plotTitle "CA1a_H2AZ_H_r1 ${region}"

    computeMatrix reference-point \
        --regionsFileName ${regionsDir}/${region}_cluster1.bed ${regionsDir}/${region}_cluster2.bed\
        --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r2_R1.bw\
        --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
        --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.bed\
        --referencePoint TSS\
        --beforeRegionStartLength 1000\
        --afterRegionStartLength 1000\
        --binSize 10\
        --missingDataAsZero\
        --sortRegions keep

    plotHeatmap \
        --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
        --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r2_R1.pdf\
        --plotTitle "CA1a_H2AZ_H_r2 ${region}"

done

export region="Cala.ShZ_vs_Cala.WT"

if [ ! -d ${deeptoolsDir}/plotHeatmap/TSS/${region}/ ]; then
    mkdir -p ${deeptoolsDir}/plotHeatmap/TSS/${region}/
fi

if [ ! -d ${deeptoolsDir}/computeMatrix/TSS/${region}/ ]; then
    mkdir -p ${deeptoolsDir}/computeMatrix/TSS/${region}/
fi

computeMatrix reference-point \
    --regionsFileName ${regionsDir}/${region}_cluster1.bed ${regionsDir}/${region}_cluster2.bed ${regionsDir}/${region}_cluster3.bed\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r1_R1.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.bed\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    --sortRegions keep

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r1_R1.pdf\
    --plotTitle "CA1a_H2AZ_H_r1 ${region}"

computeMatrix reference-point \
    --regionsFileName ${regionsDir}/${region}_cluster1.bed ${regionsDir}/${region}_cluster2.bed ${regionsDir}/${region}_cluster3.bed\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r2_R1.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.bed\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    --sortRegions keep

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r2_R1.pdf\
    --plotTitle "CA1a_H2AZ_H_r2 ${region}"

export region="LZ_vs_SZ"

if [ ! -d ${deeptoolsDir}/plotHeatmap/TSS/${region}/ ]; then
    mkdir -p ${deeptoolsDir}/plotHeatmap/TSS/${region}/
fi

if [ ! -d ${deeptoolsDir}/computeMatrix/TSS/${region}/ ]; then
    mkdir -p ${deeptoolsDir}/computeMatrix/TSS/${region}/
fi

computeMatrix reference-point \
    --regionsFileName ${regionsDir}/${region}_cluster1.bed\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r1_R1.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.bed\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    --sortRegions keep

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r1_R1.matrix.gz\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r1_R1.pdf\
    --plotTitle "CA1a_H2AZ_H_r1 ${region}"

computeMatrix reference-point \
    --regionsFileName ${regionsDir}/${region}_cluster1.bed\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_H_r2_R1.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.bed\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    --sortRegions keep

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/${region}/CA1a_H2AZ_H_r2_R1.matrix.gz\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/${region}/CA1a_H2AZ_H_r2_R1.pdf\
    --plotTitle "CA1a_H2AZ_H_r2 ${region}"


