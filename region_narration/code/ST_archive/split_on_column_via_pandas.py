#!/usr/bin/env python

import re
import os
import pandas as pd

in_fn = "../data/dhs_bed.bed"
column_index_to_split_on = 8
out_dir = "../results/split"

def main():
    df = pd.read_csv(in_fn, header=None, sep='\t')
    dfs = dict(tuple(df.groupby(df.columns[[column_index_to_split_on]].tolist())))
    skvs = sanitize_keys(dfs.keys())
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for k, v in skvs.items():
        out_fn = os.path.join(out_dir, '{}.bed'.format(v))
    dfs[k].to_csv(out_fn, index=False, sep='\t')
    
def sanitize_keys(ks):
    sks = [re.sub('[./]', '', x) for x in ks]
    sks = [re.sub('\s+', '_', x) for x in sks]
    return {y:x for (x,y) in zip(sks, ks)}

if __name__ == '__main__':
    main()


