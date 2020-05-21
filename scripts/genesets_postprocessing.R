library(deepToolsUtils)

m1 <- deepToolsUtils::computeMatrixLoader("/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/deeptools/computeMatrix/TSS/genesets_subset.matrix.gz")

m1s1means <- data.table::data.table('Gene Symbol' = m1$computeMatrixRows, 'Mean coverage (H2AZ/Input)' = rowMeans(m1$computeMatrix[,1:200]))
m1s2means <- data.table::data.table('Gene Symbol' = m1$computeMatrixRows, 'Mean coverage (H2AZ/Input)' = rowMeans(m1$computeMatrix[,201:400]))
m1s3means <- data.table::data.table('Gene Symbol' = m1$computeMatrixRows, 'Mean coverage (H2AZ/Input)' = rowMeans(m1$computeMatrix[,401:600]))

setwd('~/Data/Tremethick/Breast/Xenografts')

data.table::fwrite(m1s1means, file='MCF10Ca1a_mean_H2AZ_chip.csv')
data.table::fwrite(m1s2means, file='MCF10A_TGFb_mean_H2AZ_chip.csv')
data.table::fwrite(m1s3means, file='MCF10A_WT_mean_H2AZ_chip.csv')
