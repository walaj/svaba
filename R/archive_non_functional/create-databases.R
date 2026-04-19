####################
## load CCDS and name table
ccds <- fread("/xchip/gistic/Jeremiah/tracks/ccdsGene.hg19.txt")

## load and format name table (CCDS to gene name conversion)
nm   <- fread("/xchip/gistic/Jeremiah/tracks/name_conversion.tsv") 
setnames(nm, colnames(nm), gsub(" ","_",colnames(nm)))
app_sym <- rep(nm$Approved_Symbol, times=nm[, length(strsplit(CCDS_IDs, ",")[[1]]), by=HGNC_ID]$V1)
nmm <- data.table(ccds=gsub(" ", "", unlist(strsplit(nm$CCDS_IDs,","))), gene=app_sym)
setkey(nmm, ccds)

## merge names with CCDS
ccds[, ccds := gsub("(.*?)\\..*", "\\1", name)]
setkey(ccds, ccds)
ccds <- nmm[ccds]

## switch frames 1 and 2, so that 2 refers to 2/3 of codon, 1 to 1/3
## eg. ++***+++***+++ starts with 2/3 of codon in its sense dir, so frame is 2
ccds$exonFrames <- gsub("1","A", ccds$exonFrames)
ccds$exonFrames <- gsub("2","1", ccds$exonFrames)
ccds$exonFrames <- gsub("A","2", ccds$exonFrames)

## dedupe on gene, taking the largest one
ccds[, negc := -exonCount]
setkey(ccds, gene, negc)
ccds <- ccds[!duplicated(gene)]
ccds[, negc := NULL]

## create exon intron track
ccds[, exon_start    := min(as.numeric(strsplit(exonStarts, ",")[[1]])), by=c("name", "chrom", "txStart")]
ccds[, exon_end      := max(as.numeric(strsplit(exonEnds,  ",")[[1]])), by=c("name", "chrom", "txStart")]
ccds[, intron_ends   := { ll <- as.numeric(strsplit(exonStarts, ",")[[1]]);   paste(ll[seq(2, length(ll), length.out=ifelse(length(ll) > 1,   length(ll)-1, 0))], collapse=",")} , by=c("name", "chrom", "txStart")]
ccds[, intron_starts := { ll <- as.numeric(strsplit(exonEnds, ",")[[1]]);   paste(ll[seq(1, length(ll)-1, length.out=ifelse(length(ll) > 1, length(ll)-1, 0))], collapse=",")} , by=c("name", "chrom", "txStart")]
ccds[, intron_frame  := {
  ll <- as.numeric(strsplit(exonFrames, ",")[[1]])
  if (length(ll) == 1) ## no exons
    ""
  else if (strand == '+') ## in sense dir, give frame to exon after
    paste(ll[seq(2, length(ll))], collapse=",")
  else ## exons ordered in + direction, so sense exon after intron is actually the one BEFORE intron for this ordering
    paste(ll[seq(1, length(ll)-1)], collapse=",")
} , by=c("name", "chrom", "txStart")]
intrEnds   <- as.numeric(unlist(strsplit(ccds$intron_ends, ",")))
intrStarts <- as.numeric(unlist(strsplit(ccds$intron_starts, ",")))
intrFrame  <- as.numeric(unlist(strsplit(ccds$intron_frame, ",")))
intrName   <- rep(ccds$name, times=as.numeric(ccds$exonCount)- 1)
intrChr    <- gsub("chr", "", rep(ccds$chrom, times=as.numeric(ccds$exonCount)-1))
intrNums   <- as.numeric(ccds[, { s <- seq(length.out=ifelse(as.numeric(exonCount)==1, 0, as.numeric(exonCount)-1)); if (strand == '+') { s } else { rev(s) } }, by=c("name","chrom","txStart")]$V1)
exonNums   <- as.numeric(ccds[, { s <- seq(as.numeric(exonCount)); if (strand=='+') { s } else { rev(s) }}, by=c("name","chrom","txStart")]$V1)
intrGene   <- rep(ccds$gene, times=as.numeric(ccds$exonCount)-1)
exonGene   <- rep(ccds$gene, times=as.numeric(ccds$exonCount))
ccdsGene   <- as.numeric(ccds[, seq(as.numeric(exonCount)), by=c("name","chrom","txStart")]$V1)

ccds[, exon_end  := max(as.numeric(strsplit(exonEnds, ",")[[1]])), by=c("name", "chrom", "txStart")]
exonFrames <- as.numeric(unlist(strsplit(ccds$exonFrames, ",")))
exonStarts <- as.numeric(unlist(strsplit(ccds$exonStarts, ","))) 
exonEnds   <- as.numeric(unlist(strsplit(ccds$exonEnds, ",")))
ccdsName   <- rep(ccds$name, times=as.numeric(ccds$exonCount))
exonStrand <- rep(ccds$strand, times=as.numeric(ccds$exonCount))
#exonStarts[exonStrand == '+'] <- exonStarts[exonStrand == '+'] + 1 ## off by one?
exonStarts <- exonStarts + 1 ## off by one?
#exonEnds[exonStrand == '-'] <- exonEnds[exonStrand == '-'] + 1 ## off by one?
intrStrand <- rep(ccds$strand, times=as.numeric(ccds$exonCount)-1)
ccdschrom  <- gsub("chr","",rep(ccds$chrom, times=as.numeric(ccds$exonCount)))

gr.exons <-   sort(gr.fix(GRanges(ccdschrom, IRanges(exonStarts, exonEnds), strand=exonStrand, name=ccdsName, num=exonNums, frame=exonFrames, gene=exonGene), si))
gr.introns <- sort(gr.fix(GRanges(intrChr, IRanges(intrStarts, intrEnds), strand=intrStrand, name=intrName, num=intrNums, frame=intrFrame, gene=intrGene), si))


#### genes
gr.genes = sort(gr.fix(gr.nochr(with(fread("/xchip/gistic/Jeremiah/tracks/genes.hg19.ucsc.txt", sep="\t"), GRanges(chr, IRanges(beg, end), gene=symbol))), si))
gr.genes <- gr.genes[!grepl("^ULK|^NBPF|^MIR|^LOC|^OR|^SNO|^FAM|^ZNF|^SMN|^NF1P2|^POTEB|^RNF|^RGPD5|^RGPD2|^SNAR|^NBPF", gr.genes$gene) & width(gr.genes) < 3e6]
gr.genes$gene <- as.character(gr.genes$gene)
gg <- grbind(gr.exons, gr.introns)
#gg <- gr.val(gg, gr.genes, 'gene', sep="_")
#gg$gene[nchar(gg$gene) > 50] <- "MANYGENES"

## get complement
gg <- gr.stripstrand(gr.fix(gg, si))
gg$gene <- as.character(gg$gene)
gene.comp <- setdiff(gr.stripstrand(si2gr(si)), gg) + 1
fo <- gr2dt(gr.findoverlaps(gene.comp, gg))
fo <- fo[query.id %in% as.numeric(names(table(fo$query.id)[table(fo$query.id) == 2]))] ## now each intergenic region has front/back overlap
fo[, right := gg$gene[subject.id[2]], by=query.id]
fo[, left  := gg$gene[subject.id[1]], by=query.id]
gene.comp$left <- gene.comp$right <- ""
gene.comp$left[fo$query.id] <- fo$left
gene.comp$right[fo$query.id] <- fo$right
gene.comp <- gene.comp - 1

## tubio
ff <- fread("/xchip/gistic/Jeremiah/tracks/master_db_29062016_ranges.txt")
setnames(ff, c("V1","V2","V3","V4","V5","V6"), c("seqnames","start","end","dir","type","sample"))
ff$strand <- ifelse(ff$dir == "plus","+","-")
gr.tub <- dt2gr(ff)
