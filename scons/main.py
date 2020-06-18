import os
variables_file = 'vars.py'
vars = Variables(variables_file)

vars.Add('CIRI_SIM_EXEC', 
         'Path to the CIRI_simulator.pl file', 
         None)
vars.Add('GENE_ANNO', 
         'Input GTF formatted annotation file name', 
         '/blackhole/circrna/analyses/ccp_tuning/annotation/gencode.v29.annotation.gtf')
vars.Add('REFSEQ_DIR', 
         'Directory of reference sequence(s)', 
         '/blackhole/circrna/analyses/ccp_tuning/annotation/')
vars.Add('CIRI_SIM_OPT', 
         'Options for the CIRI_simulator.pl script', 
         None)

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
pruned_anno_cmd = 'parent_gene_fasta.py -t ${SOURCES[0]} -g ${SOURCES[1]} -b $BIOTYPE -o ' + ccp_anno_dir
pruned_anno = env.Command([os.path.join(ccp_anno_dir, f) for f in ['annotation4fasta.gtf', 
                                                                   'annotation4simulations.gtf', 
                                                                   'summary.txt']], 
                          [cirisim_trx, env['GENE_ANNO']], 
                          pruned_anno_cmd)

### the following FASTA file will be not necessary when Polyester R script 
### will consider the GTF of interest: get transcripts sequences in FASTA
#bin/gffread annotation/sim_01/lin_trx_read_genes.gtf -w annotation/sim_01/lin_trx_read_genes.fa -g /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.fa
### set FASTA in onle-line sequence entries
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < annotation/sim_01/lin_trx_read_genes.fa | tail -n +2 > annotation/sim_01/lin_trx_read_genes_oneline.fa
#
#########################
## 4. simulate linear reads from the transcripts by means of polyester
#########################
#
#./utils/simulate_linear_trx_reads.R -f annotation/sim_01/lin_trx_read_genes_oneline.fa -o reads/sim_01
#
#########################
## 5. Convert FASTA into FASTQ
#########################
#
#./utils/fasta2fastq.py reads/sim_01/sample_01_1.fasta | gzip -c > reads/sim_01/sample_01_1.fq.gz
#./utils/fasta2fastq.py reads/sim_01/sample_01_2.fasta | gzip -c > reads/sim_01/sample_01_2.fq.gz
#
#########################
## 6. concatenate circular and linear reads
#########################
#cat reads/sim_01/cirias_1.fq.gz reads/sim_01/sample_01_1.fq.gz > reads/sim_01/sim_ribo_01_1.fq.gz
#cat reads/sim_01/cirias_2.fq.gz reads/sim_01/sample_01_2.fq.gz > reads/sim_01/sim_ribo_01_2.fq.gz
#
#
#########################
## 7. prune original annotation to simulate unknown transcripts
#########################
#
### select 
### (i)   90% of genes that was used to generate circRNA reads by circ_simulator, plus 
### (ii)   5% of not expressed linear transcripts/genes, and exclude 
### (iii)  2.5% annotation that will mimic linear transcript expressed but not annotated. Also, exclude
### (iv)   2.5% annotation that will mimic expressed circRNAs with no annotation and no lin trx expressed
#
## (ii) 528+29+15=572
#grep gene_id  reads/sim_01/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -572 | tail -15 > annotation/sim_01/known_not_xprd_genes.txt
#
#cat annotation/sim_01/anno_n_xprd_genes.txt annotation/sim_01/known_not_xprd_genes.txt > annotation/sim_01/genes_to_keep.txt
#
### generate annotation file without the annotation of the genes selected above 
#grep -f annotation/sim_01/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -f annotation/sim_01/circTrx.txt > annotation/sim_01/pruned.chr1.gencode.v29.annotation.trx.gtf
#
#grep -f annotation/sim_01/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -w gene > annotation/sim_01/pruned.chr1.gencode.v29.annotation.genes.gtf
#
#cat annotation/sim_01/pruned.chr1.gencode.v29.annotation.genes.gtf annotation/sim_01/pruned.chr1.gencode.v29.annotation.trx.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/sim_01/pruned.chr1.gencode.v29.annotation.gtf 
#
#
