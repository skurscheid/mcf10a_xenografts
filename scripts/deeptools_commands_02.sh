#!/bin/bash

conda activate deeptools

export bigwigDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/bamCoverage"
export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/duplicates_removed"

# log2
bamCompare --bamfile1 ${bamDir}/A_H2AZ_H_r1_R1.DeDup.sorted.fastq_q20.bam\
           --bamfile2 ${bamDir}/MCF10A_Input_H.Q20.DeDup.sorted.bam\
           --binSize 10\
           --numberOfProcessors 32\
           --verbose\
           --outFileName MCF10A_H2AZ_log2_r1.bw

# log2
bigwigCompare --bigwig1 ${bigwigDir}/MCF10A_H2AZ_H.bw\
              --bigwig2 ${bigwigDir}/MCF10A_Input_H.bw\
              --binSize 10\
              --numberOfProcessors 32\
              --verbose\
              --outFileName MCF10A_H2AZ_over_Input.bw

# subtract
bigwigCompare --bigwig1 ${bigwigDir}/MCF10A_H2AZ_H.bw\
              --bigwig2 ${bigwigDir}/MCF10A_Input_H.bw\
              --binSize 10\
              --numberOfProcessors 32\
              --operation subtract\
              --verbose\
              --outFileName MCF10A_H2AZ_minus_Input.bw


bigwigCompare --bigwig1 ${bigwigDir}/MCF10A_TGFb_H2AZ_H.bw\
              --bigwig2 ${bigwigDir}/MCF10A_TGFb_Input_H.bw\
              --binSize 10\
              --numberOfProcessors 32\
              --verbose\
              --outFileName MCF10A_TGFb_H2AZ_over_Input.bw

bigwigCompare --bigwig1 ${bigwigDir}/MCF10CA1a_H2AZ_H.bw\
              --bigwig2 ${bigwigDir}/MCF10CA1a_Input_H.bw\
              --skipNonCoveredRegions\
              --numberOfProcessors 32\
              --verbose\
              --outFileName MCF10CA1a_H2AZ_over_Input.bw


computeMatrix reference-point --regionsFileName LZ_vs_WT.bed\
                              --scoreFileName MCF10A_H2AZ_over_Input.bw\
                              --outFileName MCF10A_H2AZ_over_Input.matrix.gz\
                              --referencePoint TSS\
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000\
                              --binSize 10\
                              --sortRegions keep

computeMatrix reference-point --regionsFileName LZ_vs_WT.bed\
                              --scoreFileName MCF10A_H2AZ_minus_Input.bw\
                              --outFileName MCF10A_H2AZ_minus_Input.matrix.gz\
                              --referencePoint TSS\
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000\
                              --binSize 10\
                              --sortRegions keep

computeMatrix reference-point --regionsFileName LZ_vs_WT.bed\
                              --scoreFileName ${bigwigDir}/MCF10A_H2AZ_H.bw\
                              --outFileName MCF10A_H2AZ.matrix.gz\
                              --referencePoint TSS\
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000\
                              --binSize 10\
                              --sortRegions keep

plotHeatmap --matrixFile MCF10A_H2AZ.matrix.gz\
            --outFileName MCF10A_H2AZ.pdf\
            --regionsLabel UP DOWN
                                                                                                             
plotHeatmap --matrixFile MCF10A_H2AZ_over_Input.matrix.gz\
            --outFileName MCF10A_H2AZ_over_Input.pdf
            
plotHeatmap --matrixFile MCF10A_H2AZ_minus_Input.matrix.gz\
            --outFileName MCF10A_H2AZ_minus_Input.pdf
                                       