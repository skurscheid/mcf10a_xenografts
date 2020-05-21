#!/bin/bash

conda activate deeptools

export bamDir="/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/duplicates_removed"

# log2
bamCompare --bamfile1 ${bamDir}/A_H2AZ_H_r1_R1.DeDup.sorted.fastq_q20.bam\
           --bamfile2 ${bamDir}/MCF10A_Input_H.Q20.DeDup.sorted.bam\
           --binSize 10\
           --numberOfProcessors 32\
           --verbose\
           --outFileName MCF10A_H2AZ_log2_r1.bw

bamCompare --bamfile1 ${bamDir}/A_H2AZ_H_r2_R1.DeDup.sorted.fastq_q20.bam\
           --bamfile2 ${bamDir}/MCF10A_Input_H.Q20.DeDup.sorted.bam\
           --binSize 10\
           --numberOfProcessors 32\
           --verbose\
           --outFileName MCF10A_H2AZ_log2_r2.bw

computeMatrix reference-point --regionsFileName LZ_vs_WT.bed\
                              --scoreFileName MCF10A_H2AZ_log2_r1.bw\
                              --outFileName MCF10A_H2AZ_log2_r1.matrix.gz\
                              --referencePoint TSS\
                              --beforeRegionStartLength 1000\
                              --afterRegionStartLength 1000\
                              --binSize 10\
                              --sortRegions keep\
                              --outFileSortedRegions sort.bed

plotHeatmap --matrixFile MCF10A_H2AZ_log2_r1.matrix.gz\
            --outFileName MCF10A_H2AZ_log2_r1.pdf\
            --regionsLabel UP DOWN
                    