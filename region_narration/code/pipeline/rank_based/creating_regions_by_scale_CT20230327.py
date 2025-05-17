import numpy as np
import sys

def main(scale, bw, chrom_len_path, out_path):
    # load lengths of chromosomes
    chrm_lengths = []
    chrms = []
    f = open(chrom_len_path, 'r')
    for line in f.readlines():
        line = line.strip().split('\t')
        chrms.append(line[0])
        chrm_lengths.append(int(line[1]))
    f.close()

    # create regions
    bed_out = []
    for c, l in zip(chrms, chrm_lengths):
        start = 0
        end = scale
        seqname = 'chr' + c
        while start < l:
            bed_out.append([seqname, start, end])
            start += bw
            end += bw
    bed_out = np.array(bed_out)

    # save
    np.savetxt(out_path, bed_out, fmt='%s', delimiter='\t')

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4])
