#!/usr/bin/env python

from __future__ import print_function
import argparse, math, re, os
from collections import defaultdict

description = """Partition input genes/transcripts and prepare annotation files

-------------------------------------------------------------------------------
case |   circRNA   | linear trx | linear trx | notes
     |  expression | expression | annotation |
---------------------------------------------|---------------------------------
  1  |  yes        | yes        | yes        | circRNA from expressed known trx
  2  |  yes        | yes        | no         | circRNA from novel trx
  3  |  yes        | no         | yes        | only circRNA from known trx
  4  |  yes        | no         | no         | novel circRNA gene
  5  |  no         | yes        | yes        | just no circRNA
  6  |  no         | yes        | no         | only novel trx
  7  |  no         | no         | yes        | known trx not expressed
  8  |  no         | no         | no         | nothing expressed, nothing known 

-- Output --
Two annotation sets/files and one summary file will be generated:
(1) annotation4fasta.gtf: annotation to get transcript FASTA that will contain cases 
    that consider linear expression (i.e. cases 1, 2, 5 and 6)
(2) annotation4simulation.gtf: annotation file to use as input in simulations that 
    will contain cases that consider known linear trx annotation, no matter whether 
    the transcript is expressed (i.e. cases 1, 3, 5 and 7)
(3) summary.txt: reporte case sets percentages and sizes, followed by each gene and 
    its relative case that has been recorded in the annotation files"""

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = description,
                                    formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--trx-list', dest= 'trxlist', required = True, type = str,
                         help = 'Gene/transcript list file')
    
    parser.add_argument('-g', '--gtf', required = True, type = str,
                         help = 'Gene annotation GTF file')

    parser.add_argument('-p', '--partition-size', dest = 'partsizes', required = False, type = str,
                         default = "85,5,5,2.5,1,1,0.5",
                         help = 'Size, in percentage, of each case set (default: "85,5,5,2.5,1,1,0.5")')

    parser.add_argument('-b', '--biotype', dest = 'gtf_tag', required = False, type = str,
                         default = 'gene', choices = ('gene', 'trx'),
                         help = 'Whether in the TRXLIST file there are gene_ids or transcript_ids (default: gene)')

    parser.add_argument('-o', '--outdir', required = False, default = 'ccp_sim',
                        type = str, 
                        help = 'Output directory (default: ccp_sim)')

    args = parser.parse_args()

    # get trx/gene list: total trx/genes?
    genes = defaultdict(int)
    with open(args.trxlist, 'r') as gene_file:
        for g in gene_file:
            genes[g.strip()] += 1

    tot_genes = len(genes.keys())

    # compute the partitions' size according to the number of input genes and input percentages
    parts = [float(p) for p in args.partsizes.split(',')]
    
    ## fill missing sizes by even quotes 
    if len(parts) < 7:
        filling_parts = 7 - len(parts)
        filling_parts_val = (100 - sum(parts)) / filling_parts
        for i in range(filling_parts):
            parts.append(filling_parts_val)

    ## compute number of genes for each case
    tot_circ_part = sum(parts[0:4])
    required_genes = int(math.ceil(tot_genes * 100 / tot_circ_part))
    part_set_sizes = [int(math.ceil(required_genes * p / 100)) for p in parts]

    ## fix part sizes
    fixed_part_set_sizes = part_set_sizes[0:3]
    if sum(part_set_sizes[0:4]) > tot_genes:
        fixed_part_set_sizes.append(part_set_sizes[3] - (sum(part_set_sizes[0:4]) - tot_genes))
    else:
        fixed_part_set_sizes.append(part_set_sizes[3])
    fixed_part_set_sizes.extend(part_set_sizes[4:7])

    ########################################################
    ## assign each gene to a case:
    ## first, assign circRNA genes to cases 1-4
    ## then, while parsing the annotation file,
    ## assign and output accordingly the non-circRNA genes
    ########################################################

    ## circRNA genes
    cumulative_sizes = [f for f in fixed_part_set_sizes]
    for i in range(1, len(cumulative_sizes)):
        cumulative_sizes[i] += cumulative_sizes[i-1]
    
    index = 0
    case = 0
    for g in genes.keys():
        index += 1
        if index > cumulative_sizes[case]:
            case += 1
        genes[g] = case + 1

    ## scan and filter annotation file according to the 
    ## partition sizes. Genes in cases 1-4 are already assigned,
    ## now will assign (and filtered along) genes in cases 5-7.
    ## There are two annotation sets/files to be generated:
    ## (1) annotation to get transcript FASTA -> will contain cases that consider linear 
    ##     expression (i.e. 1,2,5,6)
    ## (2) annotation file to use as input in simulations -> will contain cases that consider 
    ##     known linear trx annotation, no matter whether the trx is expressed (i.e. 1,3,5,7)
    
    ## keep track of how many genes is output for cases 5-7,
    ## which were not included in the genes dict
    ## start from case 5, then update counter and increase index case when counter > case size
    counter = 0
    index = 4 

    gtf_tag = 'gene_id' if args.gtf_tag == 'gene' else 'transcript_id'

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    fasta_ann = open(os.path.join(args.outdir, 'annotation4fasta.gtf'), 'w')
    simul_ann = open(os.path.join(args.outdir, 'annotation4simulations.gtf'), 'w')

    with open(args.gtf, 'r') as ann:
        for l in ann:
            if l.startswith('#'):
                continue
            ann_field = l.split('\t')[8]
            if gtf_tag in ann_field:
                gene_id = re.sub(r'.*' + gtf_tag + ' "([^"]+)";.*', r'\1', ann_field).strip()

                ## check if is a circRNA parent gene
                ## exploit the defaultdict defaults
                if genes[gene_id] == 0:
                    counter += 1
                    if index < 7:
                        if counter > fixed_part_set_sizes[index]:
                            index += 1
                            counter = 1
                        genes[gene_id] = index + 1 if index < 7 else 0

                ## (1) annotation file for FASTA
                if genes[gene_id] in (1, 2, 5, 6):
                    print(l.strip(), file = fasta_ann)
                    #'FASTA:\t' + gene_id + '\t' + str(genes[gene_id])

                ## (2) annotation file for simulation annotation
                if genes[gene_id] in (1, 3, 5, 7):
                    print(l.strip(), file = simul_ann)
                    #'ANNOT:\t' + gene_id + '\t' + str(genes[gene_id])

    fasta_ann.close()
    simul_ann.close()
    
    ## write summary file with case sets percentages and sizes
    ## followed by each gene and its relative case
    with open(os.path.join(args.outdir, 'summary.txt'), 'w') as summary_file:
        
        print(dict(zip(range(1, 8), 
                       list(zip(fixed_part_set_sizes,
                                [str(p) + '%' for p in parts])))), 
              file = summary_file)
    
        for g,c in genes.iteritems():
            if c > 0:
                print(g + '\t' + str(c), file = summary_file)

