__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-01-24"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Rules for aligning RNA-Seq reads to reference genome 
and quantifying transcript abundance with STAR
(https://github.com/alexdobin/STAR)

For usage, include this in your workflow.
"""

rule star_align:
    conda:
        "../envs/xenografts.yaml"
    version:
        "1"
    threads:
        8 
    params:
        encodeOptions = "--outFilterType BySJout --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04"
    input:
        fq1 = "fastp/trimmed/pe/{library}/{batch}/{sample}.end1.fastq.gz",
        fq2 = "fastp/trimmed/pe/{library}/{batch}/{sample}.end2.fastq.gz",
        index = lambda wildcards: config["params"]["STAR"]["index"][wildcards["ref_index"]]["gadi"]
    output:
        permDir = directory("star/{library}/{ref_index}/{batch}/{sample}/"),
        tempDir = directory("star/temp/{library}/{ref_index}/{batch}/{sample}")
    shell:
        """
            STAR --runThreadN {threads}\
                --genomeDir {input.index}\
                --readFilesIn {input.fq1} {input.fq2}\
                --outFileNamePrefix {output.permDir}\
                --outTmpDir {output.tempDir}\
                --outReadsUnmapped Fastq\
                --outSAMtype BAM Unsorted\
                --quantMode GeneCounts\
                --readFilesCommand zcat\
                {params.encodeOptions}
        """
