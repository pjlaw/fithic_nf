#!/usr/bin/env python

import numpy as np
import pandas as pd
import cooler
import sys

cooler_filepath = sys.argv[1]
chrom = sys.argv[2]
resolution = sys.argv[3]
ofile_name = sys.argv[4]

clr = cooler.Cooler(f"{cooler_filepath}::resolutions/{resolution}")

#get the bins
bins = clr.bins().fetch(chrom)
#get the raw counts - using matrix instead of pixels since fetch will only filter to region on chrom1
count_mat = clr.matrix(join=True, as_pixels=True).fetch(chrom)
#shift to midpoints - required by fithic
count_mat["mid1"] = ((count_mat["start1"]+count_mat["end1"])/2).astype(int)
count_mat["mid2"] = ((count_mat["start2"]+count_mat["end2"])/2).astype(int)

bins["weight"].replace(np.nan, -1).to_csv(f"{ofile_name}.norm", header=False, index=False, sep="\t")
count_mat[["chrom1", "mid1", "chrom2", "mid2", "count"]].to_csv(f"{ofile_name}.contactCounts.gz", header=False, index=False, sep="\t")
