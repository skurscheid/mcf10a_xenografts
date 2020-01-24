__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-01-23"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for quantifying transcript abundance with kallisto
(https://github.com/pachterlab/kallisto)

For usage, include this in your workflow.
"""

rule kallisto_quant:
    conda:
        "../envs/xenografts.yaml"
    version:
        "1"
    threads:
        16
    params:
        cli_params = config["params"]["kallisto"]["cli_params"]
    input:
        fq1 = "fastp/trimmed/pe/{library}/{batch}/{sample}.end1.fastq.gz",
        fq2 = "fastp/trimmed/pe/{library}/{batch}/{sample}.end2.fastq.gz",
        index = config["params"]["kallisto"]["index"]["mmus_hsap_ensembl99"]["gadi"]
    output:
        directory("kallisto/{library}/{batch}/{sample}")
    shell:
        """
           kallisto quant {input.fq1} {input.fq2}\
                    {params.cli_params}\
                    --index {input.index}\
                    --threads={threads}\
                    --output {output}
        """

