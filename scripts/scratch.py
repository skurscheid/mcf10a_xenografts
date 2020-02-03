import pandas as pd
wildcards = dict()
wildcards = {"batch" : "N1902403_RD_30-210828544_eukRNASEQ",
             "library" : "RNA-Seq",
             "sample" : "KDD5-2_L4",
             "ref_index" : "mmus_ensembl99"}

import yaml

# load a yaml file into dict()
stream = open('config.yaml', 'r')
config = yaml.load(stream, Loader=yaml.SafeLoader)

