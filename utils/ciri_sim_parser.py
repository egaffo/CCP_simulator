#!/usr/bin/env python

import argparse, sys
from collections import defaultdict

if __name__ == '__main__':
    
    description = '''Parse the CIRI_simulator.pl *.out output file and gives a
                     table with simulated circRNA ids with their expression represented as
                     backspliced reads'''

    parser = argparse.ArgumentParser(description = description,
                                    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', required = True, type = str,
                         help = 'The cirias.out file as output by CIRI_simulator.pl')
    
    parser.add_argument('-o', '--output', required = False, default = "-",
                        type = str, 
                        help = 'The output table with circRNA ids and relative '\
                               'backspliced read count. Default to stdout.')

    args = parser.parse_args()


    circ_ids = defaultdict(int)

    with open(args.input, 'r') as filein:
        ## assume the input file starts with a 'chr' line
        for line in filein:
            if line.startswith('chr'):
                header_line = line.split('\t')
                circ_id = header_line[-1].replace('|', '-').strip()
            if line.startswith('**'):
                circ_ids[circ_id] += 1

    fout = sys.stdout
    if args.output != '-':
        fout = open(args.output, 'w')        

    for circ_id, bks_reads in circ_ids.iteritems():
        fout.write(circ_id + '\t' + str(bks_reads) + '\n')
       
