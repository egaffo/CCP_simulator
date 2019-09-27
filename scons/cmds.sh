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

grep -w "exon\|transcript" chr1.gencode.v29.annotation.gtf | grep -f circTrx.txt > cirisimTrx.gtf
gffread cirisimTrx.gtf -w cirisimTrx.fa -g chr1.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < cirisimTrx.fa | tail -n +2 > cirisimTrx_oneline.fa

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



