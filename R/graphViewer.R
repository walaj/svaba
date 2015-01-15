require(igraph)

#dir = "/xchip/gistic/Jeremiah/TC_RUNS/HCC1143_141125/test"
dir = "/xchip/gistic/Jeremiah/TC_RUNS/HCC1143_141125/test/test250"
file = file.path(dir, "c_8:43487254-43492254_pass_0.asqg")
vert.file = paste(file, "vert", sep=".")
reads.file = file.path(dir, 'plots', 'readsForR_som.txt')
reads <- import.snowman.reads(reads.file)

#cont = 'c_8:43487254-43492254_11'
cont = 'c_8:43487254-43492254_0'
this_reads <- reads$rname[as.character(seqnames(reads)) %in% cont]

## read the edges
ln = readLines(file, n = 1000)
nskip = sum(grepl('^VT|^HT',  ln))
tab.e <- read.delim(file, skip=nskip, strings=FALSE, header=FALSE, sep=" ")
tab.e$V1 <- substring(tab.e$V1, 4)
colnames(tab.e) <- c('seq1', 'seq2', 'overlap_start1', 'overlap_end1', 'len1',
                   'overlap_start2', 'overlap_end2', 'len2', 'orientation', 'numdiff')
edges <- data.frame(from=tab.e$seq1, to=tab.e$seq2)

## read the verts
cmd = paste("grep", '"VT"', file, ">", vert.file)
system(cmd)
tab <- read.delim(vert.file, skip=1, strings=FALSE, header=FALSE, sep="\t")
colnames(tab) <- c("V", "rname", 'seq', 'ss')
verts <- data.frame(verts=unique(c(as.character(tab$rname), as.character(edges$to), as.character(edges$from))))

g <- graph.data.frame(edges, directed=TRUE, vertices=verts)


## format the verts
V(g)$names = as.character(verts$verts)
V(g)$color = "blue"
V(g)$color[V(g)$names %in% this_reads] <- 'red'


## format the edges
E(g)$lab <- tab.e$overlap_end2 - tab.e$overlap_start2

ppdf(plot(g, vertex.color=V(g)$color, vertex.label=V(g)$names, edge.label=E(g)$lab, vertex.size=2), height=40, width=40, filename='trim.pdf')
