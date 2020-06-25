#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("polyester")
suppressPackageStartupMessages(library(polyester))
suppressPackageStartupMessages(library(Biostrings))

option_list <- list(
    make_option(c("-f", "--fasta"), action = "store", type = "character",
                help = "Transcript or reference sequences in FASTA format"),
    make_option(c("-o", "--outdir"), action = "store", type = "character",
                default = "./",
                help = "Output directory were sample_xx_x.fasta will be saved"),
    make_option(c("-c", "--coverage"), action = "store", type = "numeric",
                default = 20,
                help = "Transcript coverage. Default 20 for ~20x coverage"),
    make_option(c("-g", "--gtf"), action = "store", type = "character",
                default = NULL,
                help = paste0("Gene annotation in GTF format. If this is set, ",
                              "then the --fasta must refer to reference sequences"))
)

parser <- OptionParser(usage = "%prog -f transcripts.fa -o lin_simulated_reads",
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments=F)

## prepare result dir
results.path <- arguments$outdir
dir.create(path = results.path, recursive = T, showWarnings = F)

# fasta = readDNAStringSet(arguments$fasta)
#
# fold_changes = matrix(1,
#                       ncol = 2,
#                       nrow = length(fasta))
n.trx <- sum("transcript" == scan(file = arguments$gtf,
                                  what = "character"),
             na.rm = T)
fold_changes = matrix(1,
                      ncol = 1,
                      nrow = n.trx)
#
# # ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# # here all transcripts will have ~equal FPKM
# readspertx = round(as.numeric(arguments$coverage) * width(fasta) / 100)

# simulation call:
if(!is.null(arguments$gtf)){
    simulate_experiment(seqpath = arguments$fasta,
                        gtf = arguments$gtf,
                        # reads_per_transcript = readspertx,
                        reads_per_transcript = 300,
                        num_reps = c(1, 1), #c(1),
                        fold_changes = fold_changes,
                        # meanmodel = TRUE,
                        outdir = results.path,
                        strand_specific = F,
                        paired = T)
}else{
    simulate_experiment(fasta = arguments$fasta,
                        # reads_per_transcript = readspertx,
                        reads_per_transcript = 300,
                        num_reps = c(1), #c(1, 1),
                        fold_changes = fold_changes,
                        # meanmodel = TRUE,
                        outdir = results.path,
                        strand_specific = F,
                        paired = T)
}
