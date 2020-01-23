/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_xenografts/Snakefile all_fastp_rna_seq\
    --configfile /home/150/sxk150/breast_cancer_nucleome/config.yaml\
	--use-conda\
	-d /scratch/kv78/mcf10a-xenografts\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_xenografts/cluster.json\
        --keep-going\
	-prn 
