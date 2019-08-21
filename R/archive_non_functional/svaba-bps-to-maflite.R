#!/usr/bin/env Rscript

library(optparse)
suppressMessages(suppressWarnings(require(VariantAnnotation, quietly=TRUE)))

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input bps.txt.gz file"),
    make_option(c("-o", "--output"), type = "character", default = "graph.pdf",  help = "Output MAFLITE")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

v <- readVcf(opt$input, "hg19")
inf <- info(v)
gr <- rowData(v)

## remove non indel lines
gr <- gr[inf$EVDNC=="INDEL",]
inf <- inf[inf$EVDNC=="INDEL",]

## remove unnecesry columns
#tab <- data.frame(chr=as.character(seqnames(gr)), startr=start(gr), endr=end(gr), ref=gr$REF, alt=gr$ALT)
tab <- data.frame(chr=as.character(seqnames(gr)), startr=start(gr), end=end(gr), ref=as.character(gr$REF), alt=as.character(gr$ALT), stringsAsFactors=FALSE)
del = sapply(tab$ref, nchar)  > sapply(tab$alt, nchar)
ins = sapply(tab$ref, nchar)  < sapply(tab$alt, nchar)

tab$alt[ins] <- inf$INSERTION[ins]
tab$ref[del] <- sapply(tab$ref[del], function(x) substr(x, 2, nchar(x)))

tab$alt[del] = "-"
tab$ref[ins] = "-"

tab <- cbind(tab, inf[, c("SCTG", "SPAN", "MAPQ", "TSPLIT", "NSPLIT")])
writeLines(tab, opt$output)
