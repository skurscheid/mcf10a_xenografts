import getopt, sys
from io import StringIO
from os import listdir
from os import walk
import glob
import pandas as pd

regionsDir="/home/sebastian/Data/Tremethick/Breast/Xenografts/gene_sets/subset"

myfiles = glob.glob(regionsDir + '/*.bed')
del(myfiles[1:3])

df1 = pd.read_csv(myfiles[0], sep = '\t', header=None, index_col=None)
df1 = df1.append(pd.read_csv(myfiles[1], sep = '\t', header=None, index_col=None))
df1 = df1.append(pd.read_csv(myfiles[2], sep = '\t', header=None, index_col=None))
df1 = df1.append(pd.read_csv(myfiles[3], sep = '\t', header=None, index_col=None))
df1 = df1.append(pd.read_csv(myfiles[4], sep = '\t', header=None, index_col=None))

df1.drop_duplicates(subset=3,keep='last', inplace=True)
df1.sort_values(by=[0,1], inplace=True)
df1.to_csv('meta_geneset.bed', sep='\t', header=False, index=False)