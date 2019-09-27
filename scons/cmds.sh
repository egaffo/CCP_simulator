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
perl CIRI_simulator.pl -O reads/cirias -G gencode.v29.annotation.gtf -C 20 -LC 0 -R 1 -LR 1 -L 101 -E 1 -D annotation/ -CHR1 1 -M 250 -M2 450 -PM 0 -S 70 -S2 0 -SE 0 -PSI 10

########################
# 2. get transcripts names from which circRNAs were simulated
########################

grep "^chr1" cirias.out | cut -f 3 | sort | uniq > circTrx.txt

########################
# 3. generate the transcripts FASTA(s)
########################

## (optional) select only a subset of transcripts to simulate
## linear reads in order to mimic the cases: 
## (i)  annotated and expressed linear trx
## (ii) not annotated but expressed linear trx

# (i) 75% of 590 = 442 transcripts/genes
grep gene_id  ../tools/ciri_simulator/reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -442 > anno_n_xprd_genes.txt

# (ii) 10% of 590 = 59 transcripts/genes
grep gene_id  ../tools/ciri_simulator/reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -501 | tail -59 > unknown_but_xprd_genes.txt

cat anno_n_xprd_genes.txt unknown_but_xprd_genes.txt > lin_trx_read_genes.txt

grep -w "exon\|transcript" chr1.gencode.v29.annotation.gtf | grep -f lin_trx_read_genes.txt > lin_trx_read_genes.gtf

## the following FASTA file will be not necessary when Polyester R script 
## will consider the GTF of interest
gffread lin_trx_read_genes.gtf -w lin_trx_read_genes.fa -g chr1.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < lin_trx_read_genes.fa | tail -n +2 > lin_trx_read_genes_oneline.fa

########################
# 4. simulate linear reads from the transcripts by means of polyester
########################

Rscript R_polyester/simulate_reads.R

########################
# 5. Convert FASTA into FASTQ
########################

../../../tools/fasta2fastq/fasta2fastq.py ../simulated_reads/sample_01_1.fasta | gzip -c > sample_01_1.fq.gz
../../../tools/fasta2fastq/fasta2fastq.py ../simulated_reads/sample_01_2.fasta | gzip -c > sample_01_2.fq.gz

########################
# 6. concatenate circular and linear reads
########################
cat ../ciri_simulator/cirias_1.fq.gz ../polyester/fqreads/sample_01_1.fq.gz > sim_ribo_01_1.fq.gz
cat ../ciri_simulator/cirias_2.fq.gz ../polyester/fqreads/sample_01_2.fq.gz > sim_ribo_01_2.fq.gz


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
grep gene_id  ../tools/ciri_simulator/reads/cirias.out | sed -r 's/.*gene_id "([^"]+)".*/\1/' | head -560 | tail -59 > known_not_xprd_genes.txt

cat anno_n_xprd_genes.txt known_not_xprd_genes.txt > genes_to_keep.txt

## generate annotation file without the annotation of the genes selected above 
grep -f genes_to_keep.txt chr1.gencode.v29.annotation.gtf > pruned.chr1.gencode.v29.annotation.gtf


