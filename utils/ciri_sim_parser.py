#!/usr/bin/env python

import argparse, sys
from collections import defaultdict

## ciri_simulator .out file explained: example follows.
##
## chr1    gene_id "ENSG00000203697.7"     ENST00000467384.2       chr1:223814691|223830535
## isoform1_573    223830385:223830535!-,223816364:223816482!-,223815711:223815844!-,223814691:223814859!-,
## isoform2_454    223830385:223830535!-,223815711:223815844!-,223814691:223814859!-,
## >       1       675
## >       2       676
## >       3       677
## >       4       678
## >       5       679
## **      2
## >       6       680
## **      1
## >       7       681
## >       8       682
## chr1:223814691|223830535 means the location back-spliced junction.
## 
## isoform1_573    223830385:223830535!-,223816364:223816482!-,223815711:223815844!-,223814691:223814859!-,
## 
## There are two isoforms inside this BSJ. The length of this isoform is 573nt. Exons position is listed behind and splited by ","
## 
## "-" represents for the strain.
## 
## "> 2 676" means that paired-end reads named "simulate:676" is the second pair read 
## that generated for this circRNA. but this paired-end reads do not cover the BSJ site. 
## 
## >       6       680
## **      1
## 
## "**" means that read 1 of this paired-end read can cover the BSJ site.

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

    #circ_ids = defaultdict(int)
    circ_ids = defaultdict(lambda: defaultdict(int))


    with open(args.input, 'r') as filein:
        ## assume the input file starts with a 'chr' line
        for line in filein:
            if line.startswith('chr'):
                header_line = line.split('\t')
                circ_id = header_line[-1].replace('|', '-').strip()
                circ_ids[circ_id]['gene_id'] = header_line[1].split(' ')[1].replace('"', '').strip() 
                circ_ids[circ_id]['trx_id']  = header_line[2].strip()                
            if line.startswith('**'):
                circ_ids[circ_id]['count'] += 1

    fout = sys.stdout
    if args.output != '-':
        fout = open(args.output, 'w')        

    for circ_id in circ_ids.keys():
        fout.write('\t'.join([circ_id, 
                              str(circ_ids[circ_id]['count']), 
                              circ_ids[circ_id]['gene_id'], 
                              circ_ids[circ_id]['trx_id']]) + '\n')
       
