#!/usr/bin/env Rscript

library(optparse)

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input asqg file from snowman run --write-asqg or sga"),
    make_option(c("-o", "--output"), type = "character", default = "graph.pdf",  help = "Output pdf to write the graph")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

require(igraph)

vert.file = paste(opt$input, "vert", sep=".")
##reads.file = file.path(dir, 'plots', 'readsForR_som.txt')
##reads <- import.snowman.reads(reads.file)

#cont = 'c_8:43487254-43492254_11'
##cont = 'c_8:43487254-43492254_0'
##this_reads <- reads$rname[as.character(seqnames(reads)) %in% cont]

## read the edges
ln= readLines(opt$input, n = 1000)
nskip = sum(grepl('^VT|^HT',  ln))
tab.e <- read.delim(opt$input, skip=nskip, strings=FALSE, header=FALSE, sep=" ")
tab.e$V1 <- substring(tab.e$V1, 4)
colnames(tab.e) <- c('seq1', 'seq2', 'overlap_start1', 'overlap_end1', 'len1',
                   'overlap_start2', 'overlap_end2', 'len2', 'orientation', 'numdiff')
edges <- data.frame(from=tab.e$seq1, to=tab.e$seq2)

## read the verts
cmd = paste("grep", '"VT"', opt$input, ">", vert.file)
system(cmd)
tab <- read.delim(vert.file, skip=1, strings=FALSE, header=FALSE, sep="\t")
colnames(tab) <- c("V", "rname", 'seq', 'ss')
verts <- data.frame(verts=unique(c(as.character(tab$rname), as.character(edges$to), as.character(edges$from))))

g <- graph.data.frame(edges, directed=TRUE, vertices=verts)


## format the verts
V(g)$names = as.character(verts$verts)
V(g)$color = "blue"
##V(g)$color[V(g)$names %in% this_reads] <- 'red'

## format the edges
E(g)$lab <- tab.e$overlap_end2 - tab.e$overlap_start2

pdf(opt$output, height=40, width=40)
plot(g, vertex.color=V(g)$color, vertex.label=V(g)$names, edge.label=E(g)$lab, vertex.size=2)
dev.off()

