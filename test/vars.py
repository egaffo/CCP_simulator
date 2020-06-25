##############################################################
# First, download and prepare annotation and reference files #
##############################################################

#mkdir -p annotation
#cd annotation
#wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#zcat Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | sed "s/^>1/>chr1/" > chr1.fa
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#zgrep -w "^chr1" gencode.v29.annotation.gtf.gz > chr1.gencode.v29.annotation.gtf

#################################################
# Now, CCP_SIM parameters can be set as follows #
#################################################

GENE_ANNO  = 'annotation/chr1.gencode.v29.annotation.gtf'

REFSEQ_DIR = 'annotation/'

CIRI_SIM_OPT = ['-C',   20, #coverage or max coverage (when choosing -R2) for circRNAs
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
