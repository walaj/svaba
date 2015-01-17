library(optparse)
#RLIBDIR = '/cga/meyerson/home/marcin/Software/R/x86_64-unknown-linux-gnu-library/3.1/'
#GIT.HOME = '/cga/meyerson/home/marcin/DB/git'

option_list = list(
    make_option(c("-i", "--input"), type = "character", default = NULL,  help = "Input txt file from a snowman preprocess qcreport.txt")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
    stop(print_help(parseobj))

if (!file.exists(opt$input))
  stop(paste("Input file does not exist", opt$input))

print(opt)

require(ggplot2)
require(reshape2)

## read the table
con  <- file(opt$input, open = "r")

df.mapq <- df.nm <- df.isize <- df.as <- df.xp <- df.len <- df.phred <- df.clip <- data.frame()
rg <- list()

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {

  ## found a new one
  if (grepl("READGROUP", line)) {
    thisrg = gsub("READGROUP:BI:(.*)", "\\1", line)
    rg[[thisrg]] <- data.frame(readgroup = thisrg)
  } else if (grepl("total", line)) {
    rg[[thisrg]]$total = as.numeric(gsub("total,([0-9]+)", "\\1", line))
  } else if (grepl("unmap", line)) {
    rg[[thisrg]]$unmap = as.numeric(gsub("unmap,([0-9]+)", "\\1", line))
  } else if (grepl("qcfail", line)) {
    rg[[thisrg]]$qcfail = as.numeric(gsub("qcfail,([0-9]+)", "\\1", line))
  } else if (grepl("duplicate", line)) {
    rg[[thisrg]]$duplicate = as.numeric(gsub("duplicate,([0-9]+)", "\\1", line))
  } else if (grepl("supplementary", line)) {
    rg[[thisrg]]$supp = as.numeric(gsub("supplementary,([0-9]+)", "\\1", line))    
  } else if (grepl("mapq", line)) {
    df.mapq <- rbind(df.mapq, as.numeric(strsplit(gsub("mapq,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("nm", line)) {
    df.nm <- rbind(df.nm, as.numeric(strsplit(gsub("nm,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("isize", line)) {
    df.isize <- rbind(df.isize, as.numeric(strsplit(gsub("isize,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("as", line)) {
    df.as <- rbind(df.as, as.numeric(strsplit(gsub("as,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("xp", line)) {
    df.xp <- rbind(df.xp, as.numeric(strsplit(gsub("xp,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("clip", line)) {
    df.clip <- rbind(df.clip, as.numeric(strsplit(gsub("clip,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("len", line)) {
    df.len <- rbind(df.len, as.numeric(strsplit(gsub("len,([0-9]+)", "\\1", line), ",")[[1]]))
  } else if (grepl("phred", line)) {
    df.phred <- rbind(df.phred, as.numeric(strsplit(gsub("phred,([0-9]+)", "\\1", line), ",")[[1]]))
  } else {
    stop(paste("Failed to read file at line:", line))
  }
  
}

close(con)
