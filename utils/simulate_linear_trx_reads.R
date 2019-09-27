# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("polyester")
library(polyester)
library(Biostrings)

#fasta_file <- "cirisimTrx.fa"
fasta_file <- "wg_cirisimTrx_oneline.fa"
fasta = readDNAStringSet(fasta_file)

fold_changes = matrix(1,
                      ncol = 2,
                      nrow = length(fasta))

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(fasta) / 100)

# simulation call:
#simulate_experiment('cirisimTrx.fa',
simulate_experiment(fasta_file,
                    reads_per_transcript = readspertx,
                    num_reps = c(1, 1),
                    fold_changes = fold_changes,
                    #outdir = 'simulated_reads',
                    outdir = 'wg_simulated_reads',
                    strand_specific = F,
                    paired = T)
