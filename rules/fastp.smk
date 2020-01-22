__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

def get_fastq_rna_seq(wildcards):
    """ returns fastq files for given sample_id """
    return units.loc[(wildcards["sample"], wildcards["batch"]), ["fq1", "fq2"]].dropna()

rule run_fastp_pe_rna_seq:
    conda:
        "../envs/xenografts.yaml"
    version:
        "2"
    threads:
        4
    input:
        get_fastq_rna_seq
    output:
        trimmed_read1 = "fastp/trimmed/pe/{batch}/{sample}.end1.fastq.gz",
        trimmed_read2 = "fastp/trimmed/pe/{batch}/{sample}.end2.fastq.gz",
        report_html = "fastp/report/pe/{batch}/{sample}.fastp.html",
        report_json = "fastp/report/pe/{batch}/{sample}.fastp.json"
    shell:
        """
           fastp -i {input[0]} -I {input[1]} -o {output.trimmed_read1} -O {output.trimmed_read2}\
                 --html {output.report_html}\
                 --json {output.report_json}\
                 --detect_adapter_for_pe\
                 --thread {threads}
        """

