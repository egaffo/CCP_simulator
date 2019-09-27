mkdir annotation
mkdir reads

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
perl bin/CIRI_simulator.pl -O reads/cirias -G /blackhole/circrna/analyses/ccp_tuning/annotation/gencode.v29.annotation.gtf -C 20 -LC 0 -R 1 -LR 1 -L 101 -E 1 -D /blackhole/circrna/analyses/ccp_tuning/annotation/ -CHR1 1 -M 250 -M2 450 -PM 0 -S 70 -S2 0 -SE 0 -PSI 10

## remove read mate id from read headers
sed "s_/[12]__" reads/cirias_1.fq | gzip -c > reads/cirias_1.fq.gz
sed "s_/[12]__" reads/cirias_2.fq | gzip -c > reads/cirias_2.fq.gz

########################
# 2. get transcripts names from which circRNAs were simulated
########################

grep "^chr1" reads/cirias.out | cut -f 3 | sort | uniq > annotation/circTrx.txt

########################
# 3. generate the transcripts FASTA(s)
########################

## (optional) select only a subset of transcripts to simulate
## linear reads in order to mimic the cases: 
## (i)  annotated and expressed linear trx
## (ii) not annotated but expressed linear trx

# (i) 75% of 590 = 442 transcripts/genes
grep gene_id  reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -442 > annotation/anno_n_xprd_genes.txt

# (ii) 10% of 590 = 59 transcripts/genes
grep gene_id  reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -501 | tail -59 > annotation/unknown_but_xprd_genes.txt

cat annotation/anno_n_xprd_genes.txt annotation/unknown_but_xprd_genes.txt > annotation/lin_trx_read_genes.txt

## get annotation to produce transcript fasta. NB: from full annotation select only 
## circRNA linear transcripts and discard other isoforms of the same gene to save 
## time in downstrem analysis
grep -w "exon\|transcript" /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -f annotation/lin_trx_read_genes.txt | grep -f annotation/circTrx.txt > annotation/lin_trx_read_genes.gtf

## the following FASTA file will be not necessary when Polyester R script 
## will consider the GTF of interest: get transcripts sequences in FASTA
bin/gffread annotation/lin_trx_read_genes.gtf -w annotation/lin_trx_read_genes.fa -g /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.fa
## set FASTA in onle-line sequence entries
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < annotation/lin_trx_read_genes.fa | tail -n +2 > annotation/lin_trx_read_genes_oneline.fa

########################
# 4. simulate linear reads from the transcripts by means of polyester
########################

./utils/simulate_linear_trx_reads.R -f annotation/lin_trx_read_genes_oneline.fa -o reads

########################
# 5. Convert FASTA into FASTQ
########################

./utils/fasta2fastq.py reads/sample_01_1.fasta | gzip -c > reads/sample_01_1.fq.gz
./utils/fasta2fastq.py reads/sample_01_2.fasta | gzip -c > reads/sample_01_2.fq.gz

########################
# 6. concatenate circular and linear reads
########################
cat reads/cirias_1.fq.gz reads/sample_01_1.fq.gz > reads/sim_ribo_01_1.fq.gz
cat reads/cirias_2.fq.gz reads/sample_01_2.fq.gz > reads/sim_ribo_01_2.fq.gz


########################
# 7. prune original annotation to simulate unknown transcripts
########################

## select ~50% of genes that was used to generate circRNA reads by circ_simulator
#grep gene_id  ../tools/ciri_simulator/reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -295 > genes_to_skip.txt
## generate annotation file without the annotation of the genes selected above 
#grep -v -f genes_to_skip.txt chr1.gencode.v29.annotation.gtf > pruned.chr1.gencode.v29.annotation.gtf

## select 
## (i)   75% of genes that was used to generate circRNA reads by circ_simulator, plus 
## (ii)  10% of not expressed linear transcripts/genes, and exclude 
## (iii) 10% annotation that will mimic linear transcript expressed but not annotated. Also, exclude
## (iv)   5% annotation that will mimic expressed circRNAs with no annotation and no lin trx expressed

# (ii)
grep gene_id  reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -560 | tail -59 > annotation/known_not_xprd_genes.txt

cat annotation/anno_n_xprd_genes.txt annotation/known_not_xprd_genes.txt > annotation/genes_to_keep.txt

## generate annotation file without the annotation of the genes selected above 
grep -f annotation/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -f annotation/cirisimTrx.txt > annotation/pruned.chr1.gencode.v29.annotation.trx.gtf

grep -f annotation/genes_to_keep.txt /blackhole/circrna/analyses/ccp_tuning/annotation/chr1.gencode.v29.annotation.gtf | grep -w gene > annotation/pruned.chr1.gencode.v29.annotation.genes.gtf

cat annotation/pruned.chr1.gencode.v29.annotation.genes.gtf annotation/pruned.chr1.gencode.v29.annotation.trx.gtf | sort -k1,1 -k4,4n -k5,5n > annotation/pruned.chr1.gencode.v29.annotation.gtf 


