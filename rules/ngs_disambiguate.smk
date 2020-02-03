__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-02"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Rules for disambiguating RNA-Seq reads from xenografts/mixed samples
(https://github.com/AstraZeneca-NGS/disambiguate)

For usage, include this in your workflow.
"""

def get_bams_for_disambiguate(wildcards):
    """ returns BAM files for each reference index """
    ref_index = list(config['params']['STAR']['index'].keys())
    l =[]
    for index in ref_index:
        path_list = ["samtools", wildcards["library"], index, wildcards["batch"], wildcards["sample"], "nsorted.bam"]
        l.append("/".join(path_list))
    return(l)


rule ngs_disambiguate:
    conda:
        "../envs/xenografts.yaml"
    version:
        "1"
    threads:
        2 
    params:
        cliOptions = "--aligner star"
    input:
        get_bams_for_disambiguate
    output:
        permDir = directory("ngs_disambiguate/{library}/{batch}/{sample}/")
    shell:
        """
        ngs_disambiguate {params.cliOptions}\
                         --prefix {wildcards.sample}\
                         --output-dir {output.permDir}\
                         {input[0]} {input[1]}
        """
