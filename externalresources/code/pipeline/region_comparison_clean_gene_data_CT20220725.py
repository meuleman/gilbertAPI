import os
from datetime import date
import sys

def main(genes_path, out_path):
    #cleaning up genes dataset
    f = open(genes_path, 'r')
    #finding the column names
    header = f.readline().strip().split('\t')
    #gathering data
    out_data = []
    for line in f.readlines():
        line = line.strip().split('\t')
        if (len(line) == len(header)):
            if (line[3] == ''): #when there is no gene name
                line[3] = 'nan'
            if '[' in line[4]: #removing bracket info in descriptions
                line[4] = line[4][:line[4].index('[')]
	    line[0] = 'chr' + line[0]
            line = '\t'.join(line)
            line += '\n'
            out_data.append(line)
    f.close()

    out_file = open(out_path, 'w')
    out_file.writelines(out_data)
    out_file.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
