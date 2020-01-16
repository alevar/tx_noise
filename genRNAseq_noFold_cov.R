#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(polyester,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))
suppressPackageStartupMessages(library(Biostrings,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))

# FASTA annotation
fasta_file <- args[1]
fasta = readDNAStringSet(fasta_file)

exp_file <- args[2]
exps <- as.numeric(scan(exp_file, character(), quote = ""))

readspertx = round(as.numeric(exps * width(fasta) / (as.numeric(args[3])*2)))

# simulation call:
simulate_experiment(args[1],
					fold_changes=1,
					strand_specific=TRUE,
					readlen=as.numeric(args[3]),
					paired=as.logical(as.numeric(args[5])), # args[5] should be a character anything other than 0
					error_rate=0.004,
					reads_per_transcript=readspertx,
					num_reps=c(1,0),
					outdir=args[4])