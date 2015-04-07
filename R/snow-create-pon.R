#!/usr/bin/env Rscript

library(optparse)

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = "qcreport.txt",  help = "Input file containing paths to germline "),
    make_option(c("-o", "--output"), type = "character", default = "qcreport.pdf",  help = "Output panel of normals file")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

if (!file.exists(opt$input)) {
  print(print_help(parseobj))
  stop(paste("Input file does not exist", opt$input, ". Must supply path to valid qcreport.txt file (generated from snowman preprocess or snowman run"))
}

print(opt)

## default
opt$input = '/home/unix/jwala/lung_wgs_pon_list.txt'
opt$output = "/home/unix/jwala/lung_wgs_pon.txt"

##
opt$vcf_files = read.delim(opt$input, stringsAsFactors=FALSE)[,1]
dt.all <- mclapply(opt$vcf_files, function(x) {
  print(x)
  if (!file.exists(x)) {
    warning("File does not exist")
    return (-1)
  }
  ab <- read_vcf(x)
  ab$file = x
  return(gr2dt(ab))
}, mc.cores=15)

dt <- rbindlist(dt.all)

## need to add indel tye, also svss
