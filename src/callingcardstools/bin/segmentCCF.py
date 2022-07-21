#!/usr/bin/env python

from callingCardsTools import blockify
from functools import reduce
import pandas as pd
from pybedtools import BedTool
import sys

def blocksToBED(chrom, ranges):
    output = ""
    for i in range(len(ranges) - 1):
        interval = "{}\t{}\t{}\n".format(chrom, ranges[i], ranges[i + 1])
        output += interval
    return BedTool(output, from_string = True)

if __name__ == "__main__":
    expFile = pd.read_table(sys.argv[1], header = None)
    chroms = reduce(lambda l, x: l if x in l else l+[x], expFile[0], [])
    
    for chrom in chroms:
        blocks = blockify(expFile[expFile[0] == chrom][1], algorithm="PELT")
        intervals = blocksToBED(chrom, blocks.astype(int))
        print(intervals)