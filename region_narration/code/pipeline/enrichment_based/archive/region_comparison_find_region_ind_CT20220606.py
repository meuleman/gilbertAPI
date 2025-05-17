import numpy as np
import sys
import os

def main(regions_path, genome_coords, out_path):
    #Load data
    f = open(regions_path, 'r')
    regions = []
    for line in f:
        regions.append(line.strip().split('\t'))
    f.close()
    regions = np.array(regions)
    
    #separate genomic region
    chrm = genome_coords.split(':')[0]
    start = int(genome_coords.split(':')[1].split('-')[0])
    end = int(genome_coords.split(':')[1].split('-')[1])

    #find region
    region_ind = np.where((regions[:, 0]==chrm) & 
                          (regions[:, 1].astype(int)<=start) & 
                          (regions[:, 2].astype(int)>=end))[0]

    #return error if no or more than 1 region is found
    if len(region_ind)!=1:
        sys.exit('Region not found')
    else:
        f = open(out_path, 'w')
        f.write(str(region_ind[0]))
        f.close

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
