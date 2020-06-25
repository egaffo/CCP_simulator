# CCP simulator
A tool to simulate RNA-seq reads from circular and linear RNAs

## Installation

Clone the GIT repository, enter the repo and run  

```{bash}
./utils/install.sh
```
## How to run

You will need a reference genome/chromosomes and gene annotation, for instance  

```{bash}
mkdir -p annotation
cd annotation
wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
zcat Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | sed "s/^>1/>chr1/" > chr1.fa
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
zgrep -w "^chr1" gencode.v29.annotation.gtf.gz > chr1.gencode.v29.annotation.gtf
```

Then, make your project directory and set the simulation parameters in a Python file named `vars.py`, such that in the `test` directory.  

Now, you are ready to generate the simulated reads by calling  

```{bash}
/path/to/ccp_sim/ccp_sim.sh
```

## Parameters

```
GENE_ANNO: Input GTF formatted annotation file name
    default: annotation/gencode.v29.annotation.gtf

REFSEQ_DIR: Directory of reference sequence(s)
    default: annotation/

CIRI_SIM_OPT: Options for the CIRI_simulator.pl script
    default: -C 20 -LC 0 -R 1 -LR 1 -L 101 -E 1 -CHR1 1 -M 250 -M2 450 -PM 0 -S 70 -S2 0 -SE 0 -PSI 10
```
