# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_fastp_rna_seq:
    input:
        expand("fastp/trimmed/pe/{batch}/{sample}.{ext}.fq.gz",
                batch = "N1902403_RD_30-210828544_eukRNASEQ",
                sample = ["KDD2-2_L5"],
                ext = ["end1", "end2"])

include: "rules/fastp.smk"
