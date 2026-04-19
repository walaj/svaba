#!/usr/bin/env Rscript

library(optparse)

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input asqg file from snowman run --write-asqg or sga"),
    make_option(c("-c", "--inputcc"),  type = "character", default = NULL,  help = "Input cc file from snowman run --write-asqg or sga"),  
    make_option(c("-o", "--output"), type = "character", default = "graph.pdf",  help = "Output pdf to write the graph"),
    make_option(c("-d", "--height"), type = "numeric", default = 20,  help = "Height"),
    make_option(c("-w", "--width"), type = "numeric", default = 20,  help = "Width"),
    make_option(c("-m", "--mincount"), type = "numeric", default = 1,  help = "Remove nodes with fewer than m components")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

require(igraph)
require(data.table)
#opt$input ="/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/150830/tmp.graph.after.asqg"

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
tab <- read.delim(vert.file, skip=0, strings=FALSE, header=FALSE, sep="\t")
if (ncol(tab) == 4) {
  colnames(tab) <- c("V", "rname", 'seq', 'ss')
} else if (ncol(tab) == 3) {
  colnames(tab) <- c("V", "rname", 'seq')
}
verts <- data.frame(verts=unique(c(as.character(tab$rname), as.character(edges$to), as.character(edges$from))))

g <- graph.data.frame(edges, directed=FALSE, vertices=verts)

## vert lengths
vert.lens <- structure(nchar(tab$seq), names=tab$rname)

## open the cc file
if (file.info(opt$inputcc)$size > 0) {
  cc <- fread(opt$inputcc)
  cc$col <- sample(colors(), length(unique(cc$V2)))[match(cc$V2, unique(cc$V2))]
  setkey(cc, col)
  
  V(g)$color = "black"
  for (x in cc$col) {
    V(g)$color[verts$verts %in% cc[x]$V1] <- x
  }
}

## format the verts
V(g)$names = paste(as.character(verts$verts), "len:", vert.lens[as.character(verts$verts)])

##V(g)$color[V(g)$names %in% this_reads] <- 'red'

## format the edges
E(g)$lab <- tab.e$overlap_end2 - tab.e$overlap_start2

## cluster
V(g)$community <- membership(cl<-clusters(g))
good_communities <- which(cl$csize >= opt$mincount)
g <- delete.vertices(g, V(g)[!community %in% good_communities])

pdf(opt$output, height=opt$height, width=opt$width)
plot(g, vertex.color=V(g)$color, vertex.label=V(g)$names, edge.label=E(g)$lab, vertex.size=2)
dev.off()

