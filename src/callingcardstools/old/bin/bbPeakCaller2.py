#!/bin/env python
# A (different) calling card peak caller using
# Bayesian Blocks.

import argparse
from functools import reduce
import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy.stats import poisson

# numpy log10 float max
log10FloatMax = np.log10(np.finfo(np.float64).max)

# Set up argument parser
parser = argparse.ArgumentParser(description = "Call peaks from calling card data using Bayesian blocks.", usage = "%(prog)s -p PVALUE CCF BLOCKS TTAA OUTPUT")
parser.add_argument("ccf", type = str, help = "Experiment CCF file")
parser.add_argument("blocks", type = str, help = "Bayesian Block partitioning of experiment CCF file")
parser.add_argument("ttaa", type = str, help = "BED file of TTAAs in genome")
parser.add_argument("output", type = str, help = "File to write output to")
parser.add_argument("-p", "--pValueCutoff", type = float, required = True, help = "–log10-transformed p-value cutoff")
parser.add_argument("-d", "--distance", type = int, required = False, help = "Merge features closer than this distance (bp)")
parser.add_argument("-s", "--sizeFilter", type = int, required = False, help = "Blocks larger than this will be discarded")
parser.add_argument("-r", "--refine", action = "store_true", default = False)

# Fast method for getting number of lines in a file
# For BED files, much faster than calling len() on file
# From https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == "__main__":
    args = parser.parse_args()

    ccf = BedTool(args.ccf)
    blocks = BedTool(args.blocks)
    ttaa = BedTool(args.ttaa)
    scaleFactor = file_len(args.ccf) / file_len(args.ttaa)
    
    if args.refine:
        data = blocks.intersect(ccf, wa = True, wb = True)
        df = data.to_dataframe()
        df = df.iloc[:, :6]
        df = df.rename(index = str, 
        columns = {"name": "TTAA_chrom", "score": "TTAA_start", 
                   "strand": "TTAA_end"})
        df = df.sort_values(["chrom", "start", "end", "TTAA_chrom", 
                             "TTAA_start", "TTAA_end"])
        groups = df.groupby(["chrom", "start", "end"])
        first = groups.nth(0)["TTAA_start"]
        last = groups.nth(-1)["TTAA_end"]
        joined = pd.concat([first, last], axis = 1).reset_index()
        refined = joined[["chrom", "TTAA_start", "TTAA_end"]]
        blocks = BedTool.from_dataframe(refined)

    data = blocks.intersect(ccf, c = True).intersect(ttaa, c = True)
    df = data.to_dataframe()
    df = df.rename(index = str, columns = {"name": "hops", "score": "TTAAs"})
    df["exp"] = df["TTAAs"] * scaleFactor
    df["-log10pValue"] = -np.log10(poisson.sf(df["hops"] - 1, df["exp"]))
    if args.sizeFilter:
        df["size"] = df["end"] - df["start"]
        df = df[df["size"] <= args.sizeFilter]
    outbed = BedTool.from_dataframe(df[df["-log10pValue"] >= args.pValueCutoff])
    if args.distance:
        outbed = outbed.merge(d = args.distance)
    outbed.saveas(args.output)