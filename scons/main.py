import os
variables_file = 'vars.py'
vars = Variables(variables_file)

vars.Add('CIRI_SIM_EXEC', 
         'Path to the CIRI_simulator.pl file', 
         None)
vars.Add('GENE_ANNO', 
         'Input GTF formatted annotation file name', 
         'annotation/gencode.v29.annotation.gtf')
vars.Add('REFSEQ_DIR', 
         'Directory of reference sequence(s)', 
         'annotation/')
vars.Add('CIRI_SIM_OPT', 
         'Options for the CIRI_simulator.pl script', 
         None)

vars.Add('ANNOPARTS', 
         'Fractions of annotation to assign to each circular/linear transcript annotation combination\n'\
         'case |   circRNA   | linear trx | linear trx | notes\n'\
         '     |  expression | expression | annotation |\n'\
         '---------------------------------------------|---------------------------------\n'\
         '  1  |  yes        | yes        | yes        | circRNA from expressed known trx\n'\
         '  2  |  yes        | yes        | no         | circRNA from novel trx\n'\
         '  3  |  yes        | no         | yes        | only circRNA from known trx\n'\
         '  4  |  yes        | no         | no         | novel circRNA gene\n'\
         '  5  |  no         | yes        | yes        | just no circRNA\n'\
         '  6  |  no         | yes        | no         | only novel trx\n'\
         '  7  |  no         | no         | yes        | known trx not expressed', 
         '85,5,5,2.5,1,1,0.5')

env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                  variables=vars)

Help(vars.GenerateHelpText(env))
unknown = vars.UnknownVariables()
if unknown:
    print "Unknown variables", unknown.keys()
    Exit(1)

env.SetDefault(KEEP_READ_MATE_ID = False)
env.SetDefault(CIRI_SIM_OPT = ['-C',   20, #coverage or max coverage (when choosing -R2) for circRNAs
                       '-LC',   0, #coverage or max coverage (when choosing -LR 2) for linear transcripts
                       '-R',    1, #random mode for circRNAs: 1 for constant coverage; 2 for random coverage
                       '-LR',   1, #random mode for linear transcripts: 1 for constant coverage; 2 for random coverage
                       '-L',  101, #read length(/bp) of simulated reads
                       '-E',    1, #percentage of sequencing error
                       '-CHR1', 1, #if only choose chr1 to simulate sequencing reads: 1 for yes; 0 for no
                       '-M',  250, #average(mu/bp) of insert length (major normal distribution)
                       '-M2', 450, #average(mu/bp) of insert length (minor normal distribution)
                       '-PM',   0, #percentage of minor normal distribution in total distribution
                       '-S',   70, #standard deviation(sigma/bp) of insert length
                       '-S2',   0, #standard deviation(sigma/bp) of insert length
                       '-SE',   0, #whether simulate exon skipping: 1 for yes; 0 for no
                       '-PSI', 10] #percentage of splice in for skipping exon(-SE should be 1)
                       )

env.SetDefault(BIOTYPE = 'trx') ## alternative is 'gene'

#############################
# 1. simulate circRNA reads #
#############################

#######################################
# 1.0 prepare annotation files        #
#     (to be done manually by the     #  
#     user before running this script)#
#######################################

#cd annotation
#wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#zcat Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | sed "s/^>1/>chr1/" > chr1.fa
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#gunzip gencode.v29.annotation.gtf.gz
# cd ..

##########################
# 1.1 run ciri_simulator #
##########################

ciri_read_dir = 'cirisim_reads'

cirisim_tgt = [os.path.join(ciri_read_dir, f) for f in ['cirias_1.fq.gz',
                                                        'cirias_2.fq.gz',
                                                        'cirias.out']]

cirisim_src = [env['CIRI_SIM_EXEC']] + \
              [env['GENE_ANNO'],
               env['REFSEQ_DIR']]

cirisim_cmd = ['perl ${SOURCES[0]} -O ${TARGETS[0].path.split("_1.fq.gz")[0]} '\
              '-G ${SOURCES[1]} -D ${SOURCES[2]} $CIRI_SIM_OPT']

if not env['KEEP_READ_MATE_ID']:
    ## remove read mate id from read headers
    cirisim_cmd = cirisim_cmd + \
                  ['sed -i "s_/[12]__" ${".".join(TARGETS[0].path.split(".")[0:2])}',
                   'sed -i "s_/[12]__" ${".".join(TARGETS[1].path.split(".")[0:2])}']

cirisim_cmd = cirisim_cmd +\
                   ['gzip ${".".join(TARGETS[0].path.split(".")[0:2])}',
                    'gzip ${".".join(TARGETS[1].path.split(".")[0:2])}']
              
cirisim = env.Command(cirisim_tgt,
                      cirisim_src,
                      ' && '.join(cirisim_cmd))

################################################
# 1.2 generate a 'true positive' circRNA table #
#     by parsing ciri_simulator output file    #
################################################

cirisim_pred_dir = 'cirisim_pred'

cirisim_pred_tgt = [os.path.join(cirisim_pred_dir, 'cirias_tp.csv')]
cirisim_pred_src = [cirisim[2]]
cirisim_pred_cmd = 'ciri_sim_parser.py -i ${SOURCES[0]} -o ${TARGETS[0]}'
cirisim_pred = env.Command(cirisim_pred_tgt, 
                           cirisim_pred_src,
                           cirisim_pred_cmd)


###############################################################
# 2. get transcripts names from which circRNAs were simulated #
###############################################################

cirisim_trx_tgt = [os.path.join(cirisim_pred_dir, 'circ_trx.txt')]
cirisim_trx_src = [cirisim[2]]
if env['BIOTYPE'] == 'trx':
    cut_field = '3'
else:
    cut_field = '2'
cirisim_trx_cmd = 'grep "^chr" $SOURCE | cut -f ' + cut_field + ' | sort | uniq > $TARGET'

cirisim_trx = env.Command(cirisim_trx_tgt, 
                          cirisim_trx_src,
                          cirisim_trx_cmd)


########################################
# 3. generate the transcripts FASTA(s) #
########################################
ccp_anno_dir = 'ccp_anno'
pruned_anno_cmd = 'parent_gene_fasta.py -t ${SOURCES[0]} '\
                  '-g ${SOURCES[1]} -p "$ANNOPARTS" -b $BIOTYPE -o ' + ccp_anno_dir
pruned_anno = env.Command([os.path.join(ccp_anno_dir, f) for f in ['annotation4fasta.gtf', 
                                                                   'annotation4simulations.gtf', 
                                                                   'summary.txt']], 
                          [cirisim_trx, env['GENE_ANNO']], 
                          pruned_anno_cmd)

#########################
## 4. simulate linear reads from the transcripts by means of polyester
#########################

### the following FASTA file will be not necessary when Polyester R script 
### will consider the GTF of interest: get transcripts sequences in FASTA
#bin/gffread annotation/sim_01/lin_trx_read_genes.gtf -w annotation/sim_01/lin_trx_read_genes.fa -g /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.fa
### set FASTA in onle-line sequence entries
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < annotation/sim_01/lin_trx_read_genes.fa | tail -n +2 > annotation/sim_01/lin_trx_read_genes_oneline.fa
#

# fasta_seq_cmd1 = 'gffread ${SOURCES[0]} -w ${TARGETS[0]} -g ${SOURCES[1]}'
# fasta_seq_cmd2 = '''awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < ${TARGETS[0]} | tail -n +2 > ${TARGETS[1]}'''
# fasta_seq = env.Command([], 
#                         [], 
#                         [fasta_seq_cmd1, fasta_seq_cmd2])

linreads_dir = 'lin_reads'
sim_lin_cmd = 'simulate_linear_trx_reads.R -f ${SOURCES[0]} -g ${SOURCES[1]} -o ${TARGETS[0].dir}'
sim_lin = env.Command([os.path.join(linreads_dir, f) for f in ['sample_01_1.fasta', 
                                                               'sample_01_2.fasta',
                                                               'sim_counts_matrix.rda',
                                                               'sim_rep_info.txt',
                                                               'sim_tx_info.txt']], 
                       [env['REFSEQ_DIR'], pruned_anno[0]], 
                       ' && '.join([sim_lin_cmd, 
                                    'rm ' + ' '.join([os.path.join(linreads_dir, f) for f in ['sample_02_1.fasta', 
                                                                                              'sample_02_2.fasta']
                                                     ])
                                   ])
                     )

#########################
## 5. Convert FASTA into FASTQ
#########################

fasta2fastq_cmd = 'fasta2fastq.py $SOURCE | gzip -c > $TARGET'
fasta2fastq1 = env.Command(os.path.join(linreads_dir, '${SOURCES[0].filebase}.fq.gz'), #'sample_01_1.fq.gz'
                           sim_lin[0], 
                           fasta2fastq_cmd)
fasta2fastq2 = env.Command(os.path.join(linreads_dir, '${SOURCES[0].filebase}.fq.gz'),
                           sim_lin[1], 
                           fasta2fastq_cmd)
                           
#########################
## 6. concatenate circular and linear reads
#########################
# cat reads/sim_01/cirias_1.fq.gz reads/sim_01/sample_01_1.fq.gz > reads/sim_01/sim_ribo_01_1.fq.gz
# cat reads/sim_01/cirias_2.fq.gz reads/sim_01/sample_01_2.fq.gz > reads/sim_01/sim_ribo_01_2.fq.gz

## with the following commad we also substitute semicolon characters from the
## read names, since semicolons could break some script in CCP. 
## No need for an elaborate regexp since read qualities are all ! char
cat_reads_cmd = 'zcat ${SOURCES[0]} ${SOURCES[1]} | sed "s/;/|/g" | gzip -c > $TARGET'

ccp_reads_dir = 'ccp_reads'
cat_reads_1 = env.Command(os.path.join(ccp_reads_dir, 'sim_ribodplt_1.fq.gz'), 
                          [cirisim[0], fasta2fastq1], 
                          cat_reads_cmd)
            
cat_reads_2 = env.Command(os.path.join(ccp_reads_dir, 'sim_ribodplt_2.fq.gz'), 
                          [cirisim[1], fasta2fastq2], 
                          cat_reads_cmd)
