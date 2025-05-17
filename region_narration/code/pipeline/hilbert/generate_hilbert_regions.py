import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import sys
import argparse
from time import time
import math
import os

# Coords obtained from  >>> fetchChromSizes hg38 | grep -v '_*_' | grep -v chrM
CHR_SIZES_UCSC = {"chr1": 248956422,
                  "chr10": 133797422,
                  "chr11": 135086622,
                  "chr12": 133275309,
                  "chr13": 114364328,
                  "chr14": 107043718,
                  "chr15": 101991189,
                  "chr16": 90338345,
                  "chr17": 83257441,
                  "chr18": 80373285,
                  "chr19": 58617616,
                  "chr2": 242193529,
                  "chr20": 64444167,
                  "chr21": 46709983,
                  "chr22": 50818468,
                  "chr3": 198295559,
                  "chr4": 190214555,
                  "chr5": 181538259,
                  "chr6": 170805979,
                  "chr7": 159345973,
                  "chr8": 145138636,
                  "chr9": 138394717,
                  "chrX": 156040895,
                  "chrY": 57227415}

"""
order 4 hilbert region size 1048576
order 5 hilbert region size 262144
order 6 hilbert region size 65536
order 7 hilbert region size 16384
order 8 hilbert region size 4096
order 9 hilbert region size 1024
order 10 hilbert region size 256
order 11 hilbert region size 64
order 12 hilbert region size 16
order 13 hilbert region size 4
order 14 hilbert region size 1
"""
ORDER_REGION_SIZES = {4: 1048576,
                      5: 262144,
                      6: 65536,
                      7: 16384,
                      8: 4096,
                      9: 1024,
                      10: 256,
                      11: 64,
                      12: 16,
                      13: 4,
                      14: 1}

# For various window sizes calculate region coordinates spanning the genome
def generate_hilbert_regions(order, output_dir):
    chrs_list = []
    starts_list = []
    ends_list = []
    order_check = []
    for chr in CHR_SIZES_UCSC.keys():
        order_check.append(chr)
        windows_chr = list(np.arange(0, CHR_SIZES_UCSC[chr], ORDER_REGION_SIZES[order]))
        windows_chr.append(CHR_SIZES_UCSC[chr])

        starts = windows_chr[:-1]
        ends = windows_chr[1:]
        chrs = [chr for _ in range(len(starts))]

        starts_list += starts
        ends_list += ends
        chrs_list += chrs

    if np.all(np.array(order_check) != np.array(['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY'])):
        raise ValueError("Chromosome dictionary should be in lexographical order")

    windows_df = pd.DataFrame(np.array(chrs_list).reshape(-1, 1), columns=["Chromosome"])
    windows_df["Start"] = starts_list
    windows_df["End"] = ends_list

    # output_file = output_dir / f"hilbert_order_{order}_windows.bed"
    output_file = os.path.join(output_dir, f"hilbert_order_{order}_windows.bed")
    windows_df.to_csv(output_file, sep="\t", header=None, index=False)

    return output_file


if __name__ == "__main__":
    generate_hilbert_regions(int(sys.argv[1]), sys.argv[2])
