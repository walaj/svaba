#!/usr/bin/env Rscript

## source all the required packages
source.all <- function() {

  githome <- Sys.getenv('GIT_HOME')
  if (!nchar(githome))
    githome <- '/xchip/gistic/Jeremiah/GIT'

  suppressMessages(suppressWarnings(require(ff, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(VariantAnnotation, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(rtracklayer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(plyr, quietly=TRUE)))
  print('...sourced 5 packages')
  suppressMessages(suppressWarnings(require(ggplot2, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(reshape2, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(GenomicRanges, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(popbio, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)))
  print('...sourced 10 packages')
  suppressMessages(suppressWarnings(require(bitops, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(seqinr, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Rsamtools, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(ff, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(multicore, quietly=TRUE)))
  print('...sourced 15 packages')
  suppressMessages(suppressWarnings(require(Biostrings, quietly=TRUE)))  
  suppressMessages(suppressWarnings(require(rtracklayer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(lattice, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(RColorBrewer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Matrix, quietly=TRUE)))
  print('...sourced 20 packages')
  
  source(file.path(githome, "grUtils", "grUtils.R"))
  source(file.path(githome, "trackData", "trackData.R"))
  source(file.path(githome, "marcin", "R", "functions.R"))
  source(file.path(githome, "marcin", "R", "db.R"))
  source(file.path(githome, "isva", "TaigaSig", "sigUtils.R"))
  print('...sourced grUtils.R, trackData.R, marcin/functions.R, marcin/db.R, sigUtils.R')
}

source.all()

##################
## PARSE OPTIONS
##################
library(optparse)

option_list = list(
    make_option(c("-f", "--FHworkspace"),  type = "character", default = NULL,            help = "Firehose workspace to retrieve data from"),
    make_option(c("-o", "--outdir"),  type = "character", default = getwd(),         help = "Firehose workspace to retrieve data from"),
    make_option(c("-p", "--FHpairset"),    type = "character", default = NULL,            help = "Firehose pairset to retreive data from"),
    make_option(c("-a", "--FHannotation"), type = "character", default = NULL,            help = "Firehose annoation to retreive data from"),
    make_option(c("-c", "--cores"), type = "numeric", default = 1,            help = "Number of cores to use"),
    make_option(c("-v", "--VCFlist"),      type = "character", default = "/home/unix/jwala/test.vcflist.txt",            help = "File containing a list of VCFs"))

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

### test data
opt$outdir = '/xchip/gistic/Jeremiah/tmp_sig'
dir.create(opt$outdir, showWarnings=FALSE)
opt$VCFlist = '/xchip/gistic/Jeremiah/Projects/Significance/Sanger578/list.txt'

setwd(opt$outdir)

###################
###################
