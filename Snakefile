# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
#
import os
import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

units = pd.read_csv(config["units"], sep = "\t").set_index(["sample", "batch", "library"], drop=False)

report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_fastp_rna_seq:
    input:
        expand("fastp/trimmed/pe/RNA-Seq/{batch}/{sample}.{ext}.fastq.gz",
                batch = "N1902403_RD_30-210828544_eukRNASEQ",
                sample = list(units["sample"]),
                ext = ["end1", "end2"])

rule all_kallisto_quant:
    input:
        expand("kallisto/RNA-Seq/{batch}/{sample}",
                batch = "N1902403_RD_30-210828544_eukRNASEQ",
                sample = list(units["sample"]))

rule all_star_align:
    input:
        expand("star/RNA-Seq/{ref_index}/{batch}/{sample}/",
                ref_index = ["mmus_ensembl99", "hsap_ensembl99"],
                batch = "N1902403_RD_30-210828544_eukRNASEQ",
                sample = list(units["sample"]))

include: "rules/fastp.smk"
include: "rules/kallisto.smk"
include: "rules/star.smk"
