library(sleuth)
library(data.table)
library(biomaRt)
library(tximport)

# prepare annotation for hsap/mmus mixed reference
mart.hsap <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
mart.mmus <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                              dataset = "mmusculus_gene_ensembl",
                              host = 'ensembl.org')

t2g <- data.table(biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                                "external_gene_name", "ensembl_transcript_id_version"), mart = mart.hsap))

t2g <- rbind(t2g, data.table(biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                                           "external_gene_name", "ensembl_transcript_id_version"), mart = mart.mmus)))


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version,
                          ens_gene = ensembl_gene_id, 
                          ext_gene = external_gene_name)


samples <- fread("/home/sebastian/Development/workflows/mcf10a_xenografts/samples.tsv")
units <- fread("/home/sebastian/Development/workflows/mcf10a_xenografts/units_rna-seq.tsv")

list.files("/home/sebastian/mount/gadi/mcf10a-xenografts/kallisto/RNA-Seq/N1902403_RD_30-210828544_eukRNASEQ", include.dirs = T, full.names = T)
xenograft_samples <- samples[library == "RNA-Seq" & description == "xenograft"]
xenograft_samples$directory <- paste("/home/sebastian/mount/gadi/mcf10a-xenografts/kallisto/RNA-Seq/N1902403_RD_30-210828544_eukRNASEQ/", xenograft_samples$sample, sep ="")

s2c <- dplyr::select(xenograft_samples, sample, condition = conditionB) 
s2c <- dplyr::mutate(s2c, path = xenograft_samples$directory )

so <- sleuth_prep(s2c, extra_bootstrap_summary = T, target_mapping = t2g, aggregation_column="ext_gene", num_cores=16)

kt <- kallisto_table(so)

# import data with tximport
abundances <- tximport(files = paste(xenograft_samples$directory, "abundance.h5", sep = "/"), 
                       type = "kallisto", countsFromAbundance = "scaledTPM", 
                       tx2gene = t2g[,c("target_id_version", "ext_gene")])

dTabund <- as.data.table(abundances$abundance)
colnames(dTabund) <- xenograft_samples$sample
pca1 <- ade4::dudi.pca(t(dTabund))
ade4::s.class(pca1$li, xax = 1, yax = 2, fac = as.factor(xenograft_samples$conditionB))

