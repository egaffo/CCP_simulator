mkdir annotation/sim_01
mkdir reads/sim_01

######################### 
# 1. simulate circRNA reads with ciri_simulator
#########################

## prepare annotation files
#cd annotation
#wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
#zcat Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | sed "s/^>1/>chr1/" > chr1.fa
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
#gunzip gencode.v29.annotation.gtf.gz
# cd ..
## run simulation
perl bin/CIRI_simulator.pl -O reads/sim_01/cirias -G /blackhole/circrna/analyses/ccp_tuning/annotation/gencode.v29.annotation.gtf -C 20 -LC 0 -R 1 -LR 1 -L 101 -E 1 -D /blackhole/circrna/analyses/ccp_tuning/annotation/ -CHR1 1 -M 250 -M2 450 -PM 0 -S 70 -S2 0 -SE 0 -PSI 10

## remove read mate id from read headers
sed "s_/[12]__" reads/sim_01/cirias_1.fq | gzip -c > reads/sim_01/cirias_1.fq.gz
sed "s_/[12]__" reads/sim_01/cirias_2.fq | gzip -c > reads/sim_01/cirias_2.fq.gz

## generate true positive circRNA table
./utils/ciri_sim_parser.py -i reads/sim_01/cirias.out -o reads/sim_01/cirias_tp.csv

########################
# 2. get transcripts names from which circRNAs were simulated
########################

grep "^chr1" reads/sim_01/cirias.out | cut -f 3 | sort | uniq > annotation/sim_01/circTrx.txt

########################
# 3. generate the transcripts FASTA(s)
########################

## (optional) select only a subset of transcripts to simulate
## linear reads in order to mimic the cases: 
## (i)  annotated and expressed linear trx
## (ii) not annotated but expressed linear trx

# (i) 90% of 587 = 528 transcripts/genes
grep gene_id  reads/sim_01/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -528 > annotation/sim_01/anno_n_xprd_genes.txt

# (ii) 5% of 587 = 29 transcripts/genes (520+29=557)
grep gene_id  reads/sim_01/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -557 | tail -29 > annotation/sim_01/unknown_but_xprd_genes.txt

cat annotation/sim_01/anno_n_xprd_genes.txt annotation/sim_01/unknown_but_xprd_genes.txt > annotation/sim_01/lin_trx_read_genes.txt

## get annotation to produce transcript fasta. NB: from full annotation select only 
## circRNA linear transcripts and discard other isoforms of the same gene to save 
## time in downstrem analysis
grep -w "exon\|transcript" /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -f annotation/sim_01/lin_trx_read_genes.txt | grep -f annotation/sim_01/circTrx.txt > annotation/sim_01/lin_trx_read_genes.gtf

## the following FASTA file will be not necessary when Polyester R script 
## will consider the GTF of interest: get transcripts sequences in FASTA
bin/gffread annotation/sim_01/lin_trx_read_genes.gtf -w annotation/sim_01/lin_trx_read_genes.fa -g /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.fa
## set FASTA in onle-line sequence entries
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < annotation/sim_01/lin_trx_read_genes.fa | tail -n +2 > annotation/sim_01/lin_trx_read_genes_oneline.fa

########################
# 4. simulate linear reads from the transcripts by means of polyester
########################

./utils/simulate_linear_trx_reads.R -f annotation/sim_01/lin_trx_read_genes_oneline.fa -o reads/sim_01

########################
# 5. Convert FASTA into FASTQ
########################

./utils/fasta2fastq.py reads/sim_01/sample_01_1.fasta | gzip -c > reads/sim_01/sample_01_1.fq.gz
./utils/fasta2fastq.py reads/sim_01/sample_01_2.fasta | gzip -c > reads/sim_01/sample_01_2.fq.gz

########################
# 6. concatenate circular and linear reads
########################
cat reads/sim_01/cirias_1.fq.gz reads/sim_01/sample_01_1.fq.gz > reads/sim_01/sim_ribo_01_1.fq.gz
cat reads/sim_01/cirias_2.fq.gz reads/sim_01/sample_01_2.fq.gz > reads/sim_01/sim_ribo_01_2.fq.gz


########################
# 7. prune original annotation to simulate unknown transcripts
########################

## select 
## (i)   90% of genes that was used to generate circRNA reads by circ_simulator, plus 
## (ii)   5% of not expressed linear transcripts/genes, and exclude 
## (iii)  2.5% annotation that will mimic linear transcript expressed but not annotated. Also, exclude
## (iv)   2.5% annotation that will mimic expressed circRNAs with no annotation and no lin trx expressed

# (ii) 528+29+15=572
grep gene_id  reads/sim_01/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -572 | tail -15 > annotation/sim_01/known_not_xprd_genes.txt

cat annotation/sim_01/anno_n_xprd_genes.txt annotation/sim_01/known_not_xprd_genes.txt > annotation/sim_01/genes_to_keep.txt

## generate annotation file without the annotation of the genes selected above 
grep -f annotation/sim_01/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -f annotation/sim_01/circTrx.txt > annotation/sim_01/pruned.chr1.gencode.v29.annotation.trx.gtf

grep -f annotation/sim_01/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -w gene > annotation/sim_01/pruned.chr1.gencode.v29.annotation.genes.gtf

cat annotation/sim_01/pruned.chr1.gencode.v29.annotation.genes.gtf annotation/sim_01/pruned.chr1.gencode.v29.annotation.trx.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/sim_01/pruned.chr1.gencode.v29.annotation.gtf 


