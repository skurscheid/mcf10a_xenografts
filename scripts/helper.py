from pathlib import Path
import os
import yaml
import pandas as pd


def fastp_targets(units):
    """function for creating snakemake targets for executing fastp rule"""
    t = []
    for index, row in units.iterrows():
        t.append("_".join(row['sample_id'],
                          row['library'],
                          row['condition'],
                          row['digest'],
                          str(row['replicate'])))
    return(t)
