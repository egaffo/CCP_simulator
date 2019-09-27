#!/usr/bin/env python

import argparse, sys

if __name__ == '__main__':
    
    description = '''Transform FASTA files into FASTQ format by attaching '''\
                  '''mock qualities '''

    parser = argparse.ArgumentParser(description = description,
                                    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', type = str,
                         help = 'The FASTA file to be converted. '\
                                'Set - for stdin')
    
    parser.add_argument('-o', '--output', required = False, default = "-",
                        type = str, 
                        help = 'The output name for the FASTQ file. '\
                               'Default stdout')

    args = parser.parse_args()

    f = sys.stdin
    if args.fasta != '-':
        f = open(args.fasta, 'r')

    fout = sys.stdout
    if args.output != '-':
        fout = open(args.output, 'w')        


    with f as filein:
        for line in filein:
            if line.startswith('>'):
                header = '@' + line.rstrip()[1:]
                continue
            else:
                seq = line.rstrip()
                head2 = '+'
                quals = ''.join(['!'] * len(seq))
                fout.write(header   + '\n' + \
                           seq      + '\n' + \
                           head2    + '\n' + \
                           quals    + '\n')

    fout.close()

