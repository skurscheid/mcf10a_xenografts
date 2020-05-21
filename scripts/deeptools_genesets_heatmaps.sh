#!/bin/bash
export fqDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/fastq"
export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/bowtie2"
export index="/Data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly"
export deeptoolsDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/deeptools"
export rootDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap"

conda activate deeptools

# all HALLMARKS
export regionsDir="/home/sebastian/Data/Tremethick/Breast/Xenografts/gene_sets/GeneSetTest_HALLMARK_Hs_SZtumor_vs_WT"
export regions=$(find $regionsDir -name '*bed')
export region_names=$(for i in $regions; do echo $i | awk -F/ '{print $NF}' | cut -f 1 -d '.' | awk '{gsub ("HALLMARK_", ""); print $0}'; done)

computeMatrix reference-point \
    --regionsFileName $regions\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_vs_Inp_H.bw ${deeptoolsDir}/bigwigCompare/A_TGFb_H2AZ_vs_Inp_H.bw ${deeptoolsDir}/bigwigCompare/A_H2AZ_vs_Inp_H.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/genesets.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/genesets.bed\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    --smartLabels\
    -p 32

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/genesets.matrix.gz\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/genesets.pdf

done

# subset as selected by Renae
export regionsDir="/home/sebastian/Data/Tremethick/Breast/Xenografts/gene_sets/subset"

export regions=$(find $regionsDir -name '*bed')
export region_names=$(for i in $regions; do echo $i | awk -F/ '{print $NF}' | cut -f 1 -d '.' | awk '{gsub ("HALLMARK_", ""); print $0}'; done)

computeMatrix reference-point \
    --regionsFileName $regions\
    --scoreFileName ${deeptoolsDir}/bigwigCompare/CA1a_H2AZ_vs_Inp_H.bw ${deeptoolsDir}/bigwigCompare/A_TGFb_H2AZ_vs_Inp_H.bw ${deeptoolsDir}/bigwigCompare/A_H2AZ_vs_Inp_H.bw\
    --outFileName ${deeptoolsDir}/computeMatrix/TSS/genesets_subset.matrix.gz\
    --outFileSortedRegions ${deeptoolsDir}/computeMatrix/TSS/genesets_subset.bed\
    --samplesLabel "MCF10Ca1a H2A.Z" "MCF10A TGFb H2A.Z" "MCF10A H2A.Z"\
    --referencePoint TSS\
    --beforeRegionStartLength 1000\
    --afterRegionStartLength 1000\
    --binSize 10\
    --missingDataAsZero\
    -p 32

plotHeatmap \
    --matrixFile ${deeptoolsDir}/computeMatrix/TSS/genesets_subset.matrix.gz\
    --regionsLabel ${region_names}\
    --outFileName ${deeptoolsDir}/plotHeatmap/TSS/genesets_subset.pdf
