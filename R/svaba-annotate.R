#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input VCF file"),
  make_option(c("-s", "--style"),  type = "character", default = "ncbi", help = "[nbci] for 1,2,...X, or [ucsc] for chr1, chr2, ...chrX"),
  make_option(c("-g", "--genome"),  type = "character", default = NULL, help = "genome build (hg19 or hg38)"),
  make_option(c("-b", "--db"),  type = "character", default = "refseq", help = "[refseq] gene database to use (must be refseq or gencode)"),
  make_option(c("-o", "--output"), type = "character", default = "no_id",  help = "Output annotation name")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

if (!tolower(opt$style) %in% c("ncbi","ucsc"))
  stop("Must specify style as ncbi or ucsc")

if (!tolower(opt$genome) %in% c("hg38","hg19"))
  stop("Must specify --genome as hg19 or hg38")

write("...loading required packages\n", stderr())
library(roverlaps) ## can find at https://github.com/walaj/roverlaps
##library(plyr)
library(RMySQL)
library(data.table)

## setup the connect to UCSC
write("...downloading annotation tracks from UCSC\n", stderr())
assembly <- opt$genome
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu")
query <- function(...) dbGetQuery(mychannel, ...)

if (tolower(opt$db) == tolower("refseq")) {
  ## download the exons and their gene names from RefSeq (UCSC track)
  ## - reminder that refseq ids can be interpreted:
  ##                 https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.what_is_the_difference_between
  ## - (from above): Accession numbers that begin with the prefix XM_ (mRNA), XR_ (non-coding RNA), and XP_ (protein)
  ##                 are model RefSeqs produced either by NCBIs genome annotation pipeline or copied from
  ##                 computationally annotated submissions to the INSDC.NCBI
  ##                 ...subsequently curated RefSeq records (with NM_, NR_, or NP_ accession prefixes)
  genes <- suppressWarnings(data.table::as.data.table(query(paste0("SELECT name, name2, chrom, txStart, txEnd, strand, exonStarts, exonEnds, exonCount, exonFrames FROM ", assembly, ".refGene"))))
  genes[, prefix := gsub("([A-Z]+)_.*","\\1", name)]
  genes[, longest := { a=rep(FALSE, .N); a[which.max(txEnd - txStart)] <- TRUE; a }, by=name2]
  genes[, numTranscriptAnnotations := .N, by=name2]
  genes <- genes[txEnd - txStart >= 100] ## only genes longer than 100 bp
  genes <- genes[genes$longest] ## just take the longest transcriopt
  write('...defaulting to annotating only the longest transcript')
  genes[, longest := NULL]

  setnames(genes, c("name2","name"), c("geneSymbol","id"))

} else if (tolower(opt$db) == "gencode") {

  ## if want to use GENCODE instead, need to get table to convert GENCODE ids to gene names
  genes <- suppressWarnings(data.table::as.data.table(query(paste0("SELECT name, chrom, txStart, txEnd, strand, exonStarts, exonEnds, exonCount FROM ", assembly, ".knownGene"))))
  codes <- suppressWarnings(data.table::as.data.table(query(paste0("SELECT kgID, mRNA, geneSymbol, spID, refSeq FROM ", assembly, ".kgXref"))))

  ## merge the gene ids to gene names
  setnames(genes, c("name"), c("kgID"))
  setkey(genes, kgID)
  setkey(codes, kgID)
  genes <- codes[genes]
  stopifnot(!any(duplicated(genes$kgID)))
  genes[, geneSymbol := toupper(geneSymbol)]

  setnames(genes, "kgID", "id")
} else {
  stop("db option must be either GENCODE or RefSeq")
}

write("...parsing downloaded annotation data\n", stderr())

## convert to key'ed format
setnames(genes, c("chrom", "txStart","txEnd"), c("seqnames","start","end"))
if (tolower(opt$style) == "ncbi")
  genes[, seqnames := gsub("^chr", "", seqnames)]

## make the exons structure
starts <- strsplit(genes[, exonStarts], ",")
ends   <- strsplit(genes[, exonEnds], ",")
frames   <- strsplit(genes[, exonFrames], ",")
stopifnot(sum(sapply(starts, length)) == sum(sapply(ends, length)))
stopifnot(sum(sapply(starts, length)) == sum(sapply(frames, length)))
stopifnot(sum(sapply(starts, length)) == sum(genes$exonCount))
exons <- data.table(id = genes[, rep(id, exonCount)],
                    geneSymbol = genes[, rep(geneSymbol, exonCount)],
                    exonCount = genes[, rep(exonCount, exonCount)],
                    seqnames = genes[, rep(seqnames, exonCount)],
                    frames = as.numeric(unlist(frames)),
                    start = as.numeric(unlist(starts)),
                    end = as.numeric(unlist(ends)),
                    strand = genes[, rep(strand, exonCount)])
exons[, exonNum := ifelse(strand=="+", seq(.N), rev(seq(.N))), by=id]

## make an introns track
## for the frames, define the intron from as the one of the previous exon, respecting strand
## EEEEiiiiiiiEEEEiiiiiiiiEEEEEiiiiiEEE <<< intron frame is frame of RIGHT exon (iiiiEEEE)
## EEEEiiiiiiiEEEEiiiiiiiiEEEEEiiiiiEEE >>> intron from is frame of LEFT exon (EEEEiii)
ends   <- exons[exonCount > 1, start[seq(2, .N)], by=id]
starts <- exons[exonCount > 1, end[seq(1, .N-1)], by=id]
if ("exonFrames" %in% colnames(genes)) {
  frames <- exons[exonCount > 1, {
    if (strand[1]=="+")
      frames[seq(1,.N-1)]
    else
      frames[seq(2,.N)]
  }, by=id]
  stopifnot(nrow(ends) == nrow(frames))
}
stopifnot(nrow(ends) == nrow(starts))
stopifnot(all(ends$id == starts$id))
stopifnot(min(as.numeric(names(table(exons[exonCount>1, exonCount])))) > 1)

if ("exonFrames" %in% colnames(genes)) {
  introns <- data.table(id=ends$id, start=starts$V1, end=ends$V1, frame=frames$V1)
} else {
  introns <- data.table(id=ends$id, start=starts$V1, end=ends$V1)
}
short <- exons[!duplicated(id), .(id, geneSymbol, seqnames, strand, exonCount)]
setkey(introns, id)
setkey(short, id)
introns <- short[introns]
introns[, intronCount := exonCount - 1]
introns[, exonCount := NULL]
introns[, intronNum := ifelse(strand == "+", seq(.N), rev(seq(.N))), by=id]

## make the inter-genic trackk
#ro <- roverlaps(genes, genes)
#ro[,max(end), by=query.id]$V1
#setkey(genes, seqnames, start)
#starts <- genes[,end[seq(2,.N)], by=seqnames]$V1
#ends   <- genes[,start[seq(1,.N-1)], by=seqnames]$V1

## read in file
write("...reading VCF file\n", stderr())
if (grepl("gz$", opt$input)) {
  print(paste("gunzip -c", opt$input, "| grep -v ^#"))
  vcf <- data.table::fread(cmd=paste("gunzip -c", opt$input, "| grep -v ^#"), sep="\t")
} else {
  vcf <- data.table::fread(cmd=paste("grep -v ^#", opt$input), sep="\t")
}

setnames(vcf, paste0("V", seq(9)), c("seqnames","start","id","ref","alt","qual","filter","info","geno"))

## bail early if nothing found
if (nrow(vcf) == 0) {
  write("No rearrangements found\n", stderr())
  quit(status=0)
} else {
  write(paste("Found", nrow(vcf), "rearrangements\n"), stderr())
}

## grl.ix is unique ID for each rearrangement
## grl.iix is either 1 or 2 (for each side of the rearrangement)
vcf[, grl.ix  := gsub("([0-9]+):[0-9]","\\1",id)]
vcf[, grl.iix := as.numeric(gsub("[0-9]+:([0-9])","\\1",id))]

## set the strand info for BND format
vcf[grepl("SVTYPE=BND", info), strand := ifelse(grepl("^\\[", alt) | grepl("^\\]", alt), '-', '+')]
vcf[, inv := strand[1] == strand[2], by=grl.ix]
vcf[, altstrand := rev(strand), by=grl.ix]
vcf[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", alt))]
vcf[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", alt)]
vcf[, span := ifelse(seqnames==altchr, abs(start - altpos), -1)]

## annotate with gene overlaps
ro <- roverlaps::roverlaps(vcf, genes)
ro[, geneSymbol := genes$geneSymbol[subject.id]]
ro[, geneSymbol := paste(unique(geneSymbol), collapse="_"), by=query.id]
ro[, geneTranscriptCount := genes$numTranscriptAnnotations[subject.id]]
ro[, geneStrand := genes$strand[subject.id]]
vcf[ro$query.id, geneSymbol := ro$geneSymbol]
vcf[ro$query.id, geneStrand := ro$geneStrand]
vcf[, geneSymbolAlt := geneSymbol[rev(grl.iix)], by=grl.ix]
vcf[, geneStrandAlt := geneStrand[rev(grl.iix)], by=grl.ix]
vcf[ro$query.id, geneTranscriptCount := ro$geneTranscriptCount]

## annotate with genes within 20kb
PAD <- 20e3
genes20 <- data.table::copy(genes)
genes20[, start := start - PAD/2]
genes20[, end   := end + PAD/2]

ro <- roverlaps::roverlaps(vcf, genes20)
ro[, geneSymbol20kb := genes$geneSymbol[subject.id]]
ro[, geneSymbol20kb := paste(unique(geneSymbol20kb), collapse="_"), by=query.id]
vcf[ro$query.id, geneSymbol20kb := ro$geneSymbol20kb]
vcf[, geneSymbol20kbAlt := geneSymbol20kb[rev(grl.iix)], by=grl.ix]

## overlap introns
ro <- roverlaps::roverlaps(vcf, introns)
ro[, intronNum  := introns$intronNum[subject.id]]
ro[, intronGene := introns$geneSymbol[subject.id]]
ro[, intronStart := ifelse(introns$strand[subject.id] == "+", introns$start[subject.id], introns$end[subject.id])]
ro[, intronGene := paste(unique(intronGene), collapse="_"), by=query.id]
vcf[ro$query.id, intronGene := ro$intronGene]
vcf[ro$query.id, intronNum := ro$intronNum]
vcf[ro$query.id, intronStart := ro$intronStart]
vcf[, intronGeneAlt := intronGene[rev(grl.iix)], by=grl.ix]
vcf[, intronNumAlt := intronNum[rev(grl.iix)], by=grl.ix]
vcf[, intronBpsIn  := abs(intronStart - start)]
if ("frame" %in% colnames(introns)) {
  ro[, cdsFrame := introns$frame[subject.id]]
  vcf[ro$query.id, cdsFrame  := ro$cdsFrame]
}
vcf[!is.na(intronNum), bk_msg := paste(intronBpsIn, "bp_into_intron", intronNum, "of", intronGene, sep="_")]

## overlap exons
ro <- roverlaps(vcf, exons)
ro[, exonNum  := exons$exonNum[subject.id]]
ro[, exonGene := exons$geneSymbol[subject.id]]
ro[, exonStart := ifelse(exons$strand[subject.id] == "+", exons$start[subject.id], exons$end[subject.id])]
ro[, exonGene := paste(unique(exonGene), collapse="_"), by=query.id]
vcf[ro$query.id, exonGene := ro$exonGene]
vcf[ro$query.id, exonNum := ro$exonNum]
vcf[ro$query.id, exonStart := ro$exonStart]
vcf[, exonGeneAlt := exonGene[rev(grl.iix)], by=grl.ix]
vcf[, exonNumAlt := exonNum[rev(grl.iix)], by=grl.ix]
vcf[, exonBpsIn  := abs(exonStart - start)]
vcf[!is.na(exonNum), bk_msg := paste(exonBpsIn, "bp_into_exon", exonNum, "of", exonGene, sep="_")]

## annotate fusions
vcf[, fusion := {
  f <- sum(grepl("intron|exon", bk_msg)) == 2 && (geneSymbol[1] != geneSymbol[2]);
  if (is.na(f)) { FALSE } else { f }
} , by = grl.ix]
vcf[, sense  :=
      all(geneStrand != "") &&
      ((geneStrand[1] == geneStrand[2] && strand[1] != strand[2]) || (strand[1] != strand[2] && geneStrand[1] == geneStrand[2])),
    by = grl.ix]

## DECIDE WHICH PIECE IS "first" in the fusion gene (where Tx starts)
## RAR_STRAND        GENE_STRAND
## ++ or --          must be +- or -+. Starts on +
## +- or --          msut be -- or ++. Starts on side where gstrand == rarstrand
if (any(vcf$fusion)) {
vcf[vcf$fusion, gorder := {

  mk <- c(NA,NA)
  ## inversion
  if (strand[1] == strand[2]) {
    if (geneStrand[1] == '+')
      mk <- c(2,1)
    else
      mk <- c(1,2)
  }

  ## non type
  if (geneStrand[1] == geneStrand[1])
    mk <- c(1,2)
  else
    mk <- c(2,1)

  mk

}, by=grl.ix]
} else {
  vcf[, gorder := 1] ##dummy
}

## decide if in frame
if ("cdsFrame" %in% colnames(vcf)) {
  vcf[, in_frame := {
    f = (cdsFrame[gorder[2]] - cdsFrame[gorder[1]]);
    xbases <- 0; ##xbases <- nchar(INSERTION[1]) - nchar(HOMSEQ[1]);

    if (is.na(f[1]))
      FALSE
    else if (all(grepl('intron', bk_msg)) && sense[1]) ## intron to intron
      f == 0
    else
      f + xbases %% 3 == 0
  }, by=grl.ix]

  ## make the fusion messages
  vcf[, in_frame_fusion := in_frame && fusion && sense, by=grl.ix]
  vcf[vcf$in_frame_fusion, msg := paste("In-frame sense fusion from:", bk_msg[gorder[1]], "to", bk_msg[gorder[2]]), by=grl.ix]
  vcf[!vcf$in_frame_fusion & vcf$fusion & vcf$sense, msg := paste("Out-of-frame sense fusion between:", bk_msg[1], "and", bk_msg[2]), by=grl.ix]
  vcf[!vcf$sense, msg := paste("Anti-sense fusion between:", bk_msg[1], "and", bk_msg[2]), by=grl.ix]
}

vcf[, c("fusion","gorder") := NULL]

## writing to stdout
write.table(vcf[,.(seqnames, start, start, altchr, altpos, altpos, strand, altstrand, id, ref, alt, qual, filter, info, grl.ix, grl.iix, span, geneSymbol, geneStrand, geneSymbol20kb, intronGene, intronNum, bk_msg, exonGene, exonNum, in_frame, in_frame_fusion, msg)],
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t", file=stdout())

## ##llr$bk_msg[!grepl("intron|exon", llr$bk_msg)] <- ""
## write("...annotating intergenic")
## fo <- gr.findoverlaps(llr, gene.comp)
## if (nrow(fo)) {
##   fo[, left  := gene.comp$left[subject.id],  by=query.id]
##   fo[, right := gene.comp$right[subject.id], by=query.id]
##   fo[, right_start  := end(gene.comp)[subject.id], by=subject.id]
##   fo[, left_end := start(gene.comp)[subject.id], by=subject.id]
##   fo[, left_dist  := start - left_end, by=query.id]
##   fo[, right_dist := right_start - end, by=query.id]
##   fo[, mmm := paste(left_dist, "bp to right of", left, "and", right_dist, "bp to left of", right)]
##   llr$msg[!nchar(llr$msg)] <- ""
##   llr$bk_msg[fo$query.id] <- ifelse(!nchar(llr$bk_msg[fo$query.id]), fo$mmm, llr$bk_msg[fo$query.id])
## }

## ############## make the circos plot
## write("...making Circos plot")
## library(RCircos)
## data(UCSC.HG19.Human.CytoBandIdeogram);
## chr.exclude <- NULL;
## cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
## tracks.inside <- 10;
## tracks.outside <- 0;
## RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);

## ## get the gene label dat
## #genes <- genes[width(genes) < 2e6]
## #fo1 <- gr.findoverlaps(gr1+10e3, genes)
## #fo2 <- gr.findoverlaps(gr2+10e3, genes)

## ## annoying bug with seqinfo clash on 'c'
## #if (length(fo1) && length(fo2)) {
## #  fo <- c(fo1,fo2)
## #} else if (length(fo1)) {
## #  fo <- fo1
## #} else {
## #  fo <- fo2
## #}

## write("...constructing Circos plot")

## ## set the gene labels
## gene.dat <- data.frame(Chromsome=c(dt.bks$chr1, dt.bks$chr2), chromStart=c(dt.bks$pos1, dt.bks$pos2),
##                        chromEnd=c(dt.bks$pos1,dt.bks$pos2), Gene=as.character(c(dt.bks$gene1, dt.bks$gene2)),
##                        stringsAsFactors=FALSE)
## if (nrow(gene.dat))
##   gene.dat <- gene.dat[nchar(gene.dat$Gene) > 0,]
## if (nrow(gene.dat))
##   gene.dat <- gene.dat[!duplicated(gene.dat$Gene),]

## links = data.frame()
## if (length(bks))
##   links = with(dt.bks, data.frame(Chromosome=chr1, chromStart=pos1, chromEnd=pos1, Chromsome.1=chr2, chromStart.1=pos2, chromEnd.1=pos2))

## ## plot the PDF
## pdf(file=paste0(opt$output,".pdf"), height=opt$height, width=opt$width, compress=TRUE);
## RCircos.Set.Plot.Area();
## RCircos.Chromosome.Ideogram.Plot();
## if (opt$genes != 0 && nrow(gene.dat) > 0) {
##   track.num <- 1
##   RCircos.Gene.Connector.Plot(gene.dat, track.num, "in");
##   track.num <- 2;
##   name.col <- 4;
##   RCircos.Gene.Name.Plot(gene.dat, name.col, track.num, "in");
## }
## if (nrow(links) > 0) {
##   track.num = 1
##   RCircos.Link.Plot(links, track.num, by.chromosome=TRUE) ## by.chromosome is for color
## }
## dev.off()
