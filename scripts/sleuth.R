library(sleuth)
library(data.table)

samples <- fread("~/Development/workflows/mcf10a_xenografts/samples.tsv")
units <- fread("~/Development/workflows/mcf10a_xenografts/units_rna-seq.tsv")

list.files("~/mount/gadi/mcf10a-xenografts/kallisto/RNA-Seq/N1902403_RD_30-210828544_eukRNASEQ", include.dirs = T, full.names = T)
xenograft_samples <- samples[library == "RNA-Seq" & description == "xenograft"]
xenograft_samples$directory <- paste("/home/sebastian/mount/gadi/mcf10a-xenografts/kallisto/RNA-Seq/N1902403_RD_30-210828544_eukRNASEQ/", xenograft_samples$sample, sep ="")

s2c <- dplyr::select(xenograft_samples, sample, condition = description) 
s2c <- dplyr::mutate(s2c, path = xenograft_samples$directory )

so <- sleuth_prep(s2c, extra_bootstrap_summary = F)
