/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_xenografts/Snakefile ${1}\
    --configfile /home/150/sxk150/mcf10a_xenografts/config.yaml\
	--use-conda\
	-d /scratch/kv78/mcf10a-xenografts\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_xenografts/cluster.json\
        --keep-going\
	-pr ${2}
