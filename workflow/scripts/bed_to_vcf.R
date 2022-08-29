#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(bedr)

# header is chr start end
x <- read.table(args[1], header=TRUE)
bed2vcf(x, filename=args[2], zero.based=TRUE, header=NULL, fasta=args[3])
