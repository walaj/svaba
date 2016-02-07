#!/usr/bin/env Rscript

####
# load the libraries
###
.nozzleLibraries <- function() {

  print('...loading libraries')
  RLIBDIR = '/xchip/gistic/Jeremiah/R/x86_64-unknown-linux-gnu-library/3.1/'
  GIT.HOME = '/xchip/gistic/Jeremiah/GIT/'
  ISVA.HOME = paste(Sys.getenv('GIT_HOME'),  '/isva/', sep = "");
  .libPaths(c(.libPaths(), RLIBDIR))
  suppressMessages(suppressWarnings(require(methods, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(optparse, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Rsamtools, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Matrix, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(bitops, quietly=TRUE)))
    suppressMessages(suppressWarnings(require(GenomicAlignments, quietly=TRUE)))
    
  #source(file.path(GIT.HOME,"isva", "sigUtils.R"))
  #source(file.path(GIT.HOME,"isva", "sigUtils.R"))
  source(file.path(GIT.HOME,"dev/gChain/", "gChain.R"))
  ##source(file.path(GIT.HOME,"grUtilts", "grUtils.R"))
  source(file.path(GIT.HOME,"dev/gTrack/gTrack/R/", "gTrack.R"))  
  suppressMessages(suppressWarnings(require(Rsamtools, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(RColorBrewer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Nozzle.R1, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)))
  #source('/home/unix/jwala/GIT/isva/grUtils.R')
  #source('/home/unix/jwala/GIT/isva/trackData.R')

}

###
# get the base filename
###
file.name = function(paths) 
  return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths))

###
# load the breakpoint files
###
.loadBreakpoints <- function(opt) {

  ## test
  #fp <- file.path("/xchip/gistic/Jeremiah/Projects/SnowmanWithIndel_140218","bps.txt.gz")
  #tab = read.table(gzfile(fp), header=T, sep='\t', stringsAsFactors=FALSE)
  #tab$reads = ""

  gr <- GRanges(c(tab$chr1, tab$chr2), IRanges(c(tab$pos1,tab$pos2), width=1)) #, strand=ifelse(c(tab$strand1, tab$strand2)=='+', 1,-1))
  strand(gr) = c(tab$strand1, tab$strand2)
  grl <- split(gr, tab$contig)

  ## debug
  tabr = table(tab$contig)
  grl = grl[names(grl) %in% names(tabr[tabr==1])]
  mcols(grl) = tab[tab$contig %in% names(tabr[tabr==1]), ]

  return(list(som=grl, ger=GRangesList()))
  
  ## load the breakpoints
  #breakpoints.somatic.file = file.path(opt$indir, 'breakpoints', 'breakpoints.somatic.txt')
  #print(paste('Loading the breakpoints:', breakpoints.somatic.file))
  #grl.sno <- sig.load.snow(breakpoints.somatic.file)[[1]]
  
  #breakpoints.germline.file = file.path(opt$indir, 'breakpoints', 'breakpoints.germline.txt')
  #print(paste('Loading the breakpoints:', breakpoints.germline.file))
  #grl.sno.g <- sig.load.snow(breakpoints.germline.file)[[1]]
  #return(list(som=grl.sno, ger=grl.sno.g))
}

###
#
###
.initializeReport <- function(r) {

  intro.str = "Somatic Analysis -- SnowmanSV (Somatic/germline structural variation detction by string graph assembly, Wala et al in prep) is a tool for detecting germline and somatic structural variations (SVs) using whole genome sequencing (WGS) data. The engine of Snowman is a modified version of String Graph Assembler (SGA) [Simpson & Durbin, Genome Research 2012]. Snowman peforms a rolling tiled assembly across the genome in slightly overlapping 5kb windows. For efficiency, only reads which are clipped, unmapped or discordant are used in the assembly. These assembled contigs are then aligned to the genome with BWA-MEM [Heng Li, 2013], and are scored based on their alignment properties and the number of supporting sequencing reads. The final result is a list of junctions and their read support."
  if (opt$germline)
    intro.str = "Germline Analysis -- SnowmanSV (Somatic/germline structural variation detction by string graph assembly, Wala et al in prep) is a tool for detecting germline and somatic structural variations (SVs) using whole genome sequencing (WGS) data. The engine of Snowman is a modified version of String Graph Assembler (SGA) [Simpson & Durbin, Genome Research 2012]. Snowman peforms a rolling tiled assembly across the genome in slightly overlapping 5kb windows. For efficiency, only reads which are clipped, unmapped or discordant are used in the assembly. These assembled contigs are then aligned to the genome with BWA-MEM [Heng Li, 2013], and are scored based on their alignment properties and the number of supporting sequencing reads. The final result is a list of junctions and their read support."
  files.str = "Snowman produces a series of output files, which are useful for visualizing and understanding the junctions. Two core bam files are produced at the end of 'snowman run' -- r2c_clean.bam and all_bwa.bam. r2c_clean.bam is a BAM containing only reads with mappings to assembled contigs. Reads which do not assemble into contigs are discarded. Each read is in its originally aligned position, and is given a CN tag (contig) and AL tag (alignment-on-contig). all_bwa.bam is an alignment of all of the contigs to the human reference genome. r2c_final.bam and contigs_final.bam are trimmed versions of these which contain only high-confidence read-to-contig mappings, and only contigs which support either a germline or somatic event. In the 'breakpoints' folder, there are VCFs and txt files providing the breakpoints, both somatic and germline. The 'alignments' folder contains ASCII plots of the contigs and their read support. The 'plots' folder is for internal use (required to make nozzle output)."
  
  r <- setReportSubTitle(r, "Jeremiah Wala")
  r <- addToIntroduction(r, newParagraph(intro.str))
  r <- addToIntroduction(r, newParagraph(files.str))
  
  #r <- addToSummary(r, newParagraph(intro.str))
  return (r)
}

###
#
###
.overviewFigure <- function(grl.sno, opt, r) {

  print('...making the overview trackData figure')

  grsno = grl.sno$som
  suffix = 'somatic'
  if (opt$germline) {
    grsno = grl.sno$ger
    suffix = 'germline'
  }

  td <- trackData()
  td$xaxis.newline=FALSE
  td$xaxis.chronly=TRUE
  td$xaxis.nticks = 0
  overview.fig = file.path(opt$outdir, paste('breakpoint_overview_', suffix, '.png', sep=''));
  png(overview.fig, width=6000, height=800)
  display(td, links=grsno, windows=gr.all())
  suppressMessages(dev.off())
  fig = newFigure(file.name(overview.fig), fileHighRes = file.name(overview.fig), exportId = "FIGURE_1", "Track view of entire JaBbA graph, with unincorporated rearrangements removed.   Segments are assigned a copy number and plotted on the y axis.  Reference edges are shown in very light gray.  Red edges represent aberrant junctions and their weight corresponds to their inferred copy number. Loose ends are shown as dangling dotted blue edges.  The bottom track shows genome-wide normalized coverage")
  r = addToResults(r, fig)
  return (r)
  
}

###
#
###
.makeCircos <- function(grl.sno, opt, r) {

  if (!opt$germline) {
    
    ## add the somatic circos figure
    circos.fig <- file.path(opt$outdir, 'circos_somatic.png')
    circos.conf <- file.path(opt$outdir, 'circos_somatic.conf')
    circos.svg <- file.path(opt$outdir, 'circos_somatic.svg')
    gr2circos(opt$outdir, ra=grl.sno$som)
    system(paste('cd', opt$outdir, '; /home/unix/marcin/Software/circos/circos-0.62-1/bin/circos'))
    system(paste('mv', file.path(opt$outdir, 'circos.png'), circos.fig))
    system(paste('mv', file.path(opt$outdir, 'circos.svg'), circos.svg))
    system(paste('mv', file.path(opt$outdir, 'circos.conf'), circos.conf))
    fig = newFigure(file.name(circos.fig), fileHigRes=file.name(circos.fig), exportId = "FIGURE_2", "Somatic events")
    r = addToResults(r, fig)
    
  } else {

    ## add the germline circos figure
    circos.conf <- file.path(opt$outdir, 'circos_germline.conf')
    circos.svg <- file.path(opt$outdir, 'circos_germline.svg')
    circos.fig <- file.path(opt$outdir, 'circos_germline.png')
    gr2circos(opt$outdir, ra=grl.sno$ger)
    system(paste('cd', opt$outdir, '; /home/unix/marcin/Software/circos/circos-0.62-1/bin/circos'))
    system(paste('mv', file.path(opt$outdir, 'circos.png'), circos.fig))
    system(paste('mv', file.path(opt$outdir, 'circos.svg'), circos.svg))
    system(paste('mv', file.path(opt$outdir, 'circos.conf'), circos.conf))
    fig = newFigure(file.name(circos.fig), fileHigRes=file.name(circos.fig), exportId = "FIGURE_2", "Germline events")
    r = addToResults(r, fig)
  }
  
  return(r)
  
}

###
#
###
.getContigs <- function(opt) {
  
  ## read in the contigs
  snow.contigs <- file.path(opt$indir, 'contigs.bam')
  if (!file.exists(snow.contigs))
    snow.contigs <- file.path(opt$indir, 'all_bwa.bam')
  if (!file.exists(snow.contigs))
    stop('Cannot find contigs bam file')
  print(paste('...importing contigs from:', snow.contigs))
  gr.contigs <- import.snowman.contigs(snow.contigs, paste(snow.contigs, "bai", sep='.')) 

  uniq = unique(mcols(grl.sno$som)$contig)
  if (opt$germline)
    uniq = unique(mcols(grl.sno$ger)$cname)
  gr.contigs <- gr.contigs[gr.contigs$qname %in% uniq] ## keep only contigs in somatic breaks
  suppressWarnings(gr.contigs <- gr.contigs[order(nchar(gr.contigs$seq), decreasing=TRUE)]) ## put longest seqs in front, to deal with potential hard clipping
  
  ix <- !is.na(gr.contigs$cigar)
  gr.contigs <- gr.contigs[ix]
  return(gr.contigs)
  
}

###
#
###
.getReads <- function(opt) {
  
  ## read in the reads
  #snow.reads <- file.path(opt$indir, 'plots/readsForR_som.txt')
  #if (opt$germline)
  #    snow.reads <- file.path(opt$indir, 'plots/readsForR_ger.txt')
  #if (!file.exists(snow.reads))
  #  stop(paste('Cannot find', snow.reads))
  #print(paste('`...importing reads from:', snow.reads))
  #gr.reads   <- import.snowman.reads(snow.reads)

  bam = file.path(opt$indir, "r2c_clean.bam")
  bami = file.path(opt$indir, "r2c_clean.bam.bai")
  gr.reads = read.bam(bam, bami=bami, tag=c("AL", "SW", "CN", "SR", "TS"), pairs.grl=FALSE)
  #gr.reads$tn = substring(gr.reads$rname, 1,1)
  #gr.reads$rqname = gsub("[a-z]+[0-9]+_?(.*)", "\\1", gr.reads$rname)
  #gr.reads$rheader = gsub("[^_]+($)", "\\1", gr.reads$rname)
  #gr.reads$rqname = substring(gr.reads$rname, nchar(gr.reads$rheader)+1, nchar(gr.reads$rname))
  #gr.reads$flag = as.numeric(gsub("[a-z]+([0-9]+)_?(.*)", "\\1", gr.reads$rheader))
  gr.reads <- gr.reads[!is.na(gr.reads$AL) & !is.na(gr.reads$CN) & !is.na(gr.reads$SR)]
  al <- unlist(strsplit(gr.reads$AL, "x"))
  cn <- unlist(strsplit(gr.reads$CN, "x"))
  len <- sapply(strsplit(gr.reads$AL, "x"), length)
  seq <- rep(gr.reads$TS, len)
  seq[is.na(seq)] = gr.reads$seq[is.na(seq)]

  #ix = !is.na(as.integer(al)) ##debug
  grr <- GRanges(seqnames=cn, IRanges(as.integer(al), width=nchar(seq)), seq=seq)
  grr$cigar = paste(nchar(seq), "M", sep="")
  grr$rname = rep(gr.reads$qname, len)
  
  return(grr)
}

###
#
###

.getTrackData <- function(cgc, r2g, gr.contigs.this, gr.reads.this) {

  ## subet them
  cgc_t <- gSubset(cgc, xnames=unique(gr.contigs.this$qname))
  r2g_t <- gSubset(r2g, xnames=gr.reads.this$rname)
  if (length(links(r2g_t)$x)  == 0 || length(links(cgc_t)$x) == 0)
    return (trackData())

  umap <- gr.reads.this$rname[bitAnd(gr.reads.this$flag, 4) != 0]
  
  ix = !duplicated(gr.contigs.this$qname)
  cseq.set <- DNAStringSet(gr.contigs.this$seq[ix])
  names(cseq.set) <- gr.contigs.this$qname[ix]
  ix <- !duplicated(gr.reads.this$rname)
  rseq.set <- DNAStringSet(gr.reads.this$seq[ix])
  names(rseq.set) <- gr.reads.this$rname[ix]
  
  grl.contig.seq <- seq2grl(cseq.set) ## max because of hard clipping issue
  grs.contig <- lift(cgc_t, grl.contig.seq)
  td.c2g <- do.call('trackData', c(list(data=grs.contig, track.name='Contig', labels.suppress=TRUE)))
  
  grs.reads <- seq2grl(rseq.set, sn=names(rseq.set))
  reads2genome <- lift(r2g_t, grs.reads, pintersect=TRUE)
  r2g.umap <- reads2genome[ names(reads2genome) %in% umap]
  r2g.map  <- reads2genome[!names(reads2genome) %in% umap]
  td.r2g.u <- do.call('trackData', c(list(data=r2g.umap, track.name='UMap', labels.suppress=TRUE)))
  td.r2g.m <- do.call('trackData', c(list(data=r2g.map, track.name='Map', labels.suppress=TRUE)))
  td.r2g <- c(td.r2g.u, td.r2g.m)
  
   ## lift the contig onto the genome
  gr.cdum <- GRanges(names(cseq.set),IRanges(1,width(cseq.set)),strand='+')
  suppressWarnings(rr <- lift(cgc_t, gr.cdum))
  rr$mapq <- values(cgc)$mapq[rr$link.id]
  rr$col[!is.na(rr$mapq)] <- alpha('darkred', rr$mapq[!is.na(rr$mapq)]/60)
  rr$col[ is.na(rr$mapq)] <- 'yellow'
  rr$border <- 'black'
  rr.split <- split(rr, rr$query.id)
  names(rr.split) = names(cseq.set)
  
  ## make contig mapq track data
  td.c <- do.call('trackData', c(list(data=rr.split, labels.suppress=TRUE, draw.paths=TRUE)))
  
  td.r2g$height = 10
  td.c2g$height = 1
  td.c$height = 1
  
  tds <- c(td.r2g, td.c2g, td.c)
  tds$xaxis.nticks=2
  tds$xaxis.newline=FALSE
  tds$xaxis.cex = 2.0
  tds$xaxis.unit = 1;
  tds$xaxis.suffix = ""
  return(tds)
  
}

.getN50 <- function(opt, r) {

  ## load the contigs
  print('...reading contig fasta')
  contigs.file = file.path(opt$indir, 'all_contigs_bootstrap.fa')
  if (!file.exists(contigs.file)) {
    print(paste("Contigs file for estimating N50 does not exist:", contigs.file))
    return (r)
  }
    
  awk.cmd = paste("awk 'NR % 2 == 0 {print length}'", contigs.file)
  ab <- as.numeric(read.table(pipe(awk.cmd),header=FALSE)$V1)
  n50 <- Biostrings::N50(ab)

  ymax = max(table(cut(log(ab,10),breaks=seq(0,4,by=0.01))))
  lab = seq(from=2, to=3.5, by=0.1)
  p <- ggplot() + geom_histogram(data=data.frame(len=ab), aes(x=log(len,10)), binwidth=0.01, color='black', fill=NA) +
    geom_line(data=data.frame(x=rep(log(n50,10),2), y=c(0, ymax)), aes(x=x, y=y), color='red') +
    geom_text(data=data.frame(x=0.1 + log(n50,10), y=ymax), aes(x=x, y=y, label=paste("N50: ", n50, 'bp', sep='')), color='red') +
    theme_bw() + xlab('Contig length') + ylab('Count') + scale_x_continuous(breaks=lab, labels=parse(text=paste('10', lab, sep='^')))

  png(file.path(opt$outdir, 'N50hist.png'), width=800, height=400)
  print(p)
  dev.off()

  fig = newFigure('N50hist.png', fileHigRes='N50hist.png', exporId="FIGURE_1000",
    paste("Distribution of lengths of assembled contigs. The N50 is", n50))
  r = addToResults(r, fig)
  return (r)
  
}

.getTrackDataAlign <- function(cgc, r2g, gr.contigs.this, gr.reads.this) {

  cgc_t <- gSubset(cgc, xnames=unique(gr.contigs.this$qname))
  r2g_t <- gSubset(r2g, xnames=gr.reads.this$rname)
  if (length(links(r2g_t)$x)  == 0 || length(links(cgc_t)$x) == 0)
    return (trackData())

  gr.reads.this$col <- 'gray'
  gr.reads.this$border <- 'black'
  
  ## isolate discordant read pairs
  tab <- table(gr.reads.this$rname)
  gr.reads.this$col[gr.reads.this$rname %in% names(tab[tab==2])] <- 'blue'

  ## isolate the unmapped reads
  gr.reads.this$border[bitAnd(gr.reads.this$flag, 4) != 0] <- 'red'

  ## make the new reads struct
  grr <- GRanges(gr.reads.this$rname, IRanges(1, nchar(gr.reads.this$seq)))
  grr$col <- gr.reads.this$col
  grr$border <- gr.reads.this$border
  grr$rname <- gr.reads.this$rname

  grl.grr <- split(grr, grr$rname)
  reads2genome <- lift(r2g_t, grl.grr)
  
  tdo <- trackData(reads2genome, labels.suppress=TRUE)
  return(tdo)
}

#####################################
########### RUN NOZZLE ##############
#####################################

## load the libraries
.nozzleLibraries()

option_list = list(
  make_option(c("-i", "--indir"), type = "character", default = './', help = "Path to breakpoints file breakpoints.somatic.text in Snowman output directory, or root Snowman output directory"),
  make_option(c("-g", "--germline"), type = "character", default = 'false', help = "Set this flag to true to output germline snowman results instead of somatic. Default: false"),
  make_option(c("-m", "--maxpngs"), type = "integer", default = 1e8, help = "Maximum PNGS to produce per report subheading"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, help = "Number of cores to use"),
  make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

## make the output directory
if (!grepl("nozzle", opt$outdir)) 
  opt$outdir = file.path(opt$outdir, 'nozzle')
dir.create(opt$outdir, showWarnings=FALSE)

opt$germline = grepl("t|T", opt$germline)

## load the breakpoints
grl.sno = .loadBreakpoints(opt)

##debug
grl.sno$ger = grl.sno$som
# Phase 1: create report elements
r <- newReport( "SnowmanSV -- Structural variation detection by genome-wide local assembly" )
r <- .initializeReport(r)

## add the overview figure
r <- .overviewFigure(grl.sno, opt, r)

## add the contig N50 figure
#r <- .getN50(opt, r)

## make the circos figures
r <- .makeCircos(grl.sno, opt, r)

## grab the contigs and reads
gr.contigs <- .getContigs(opt)
gr.reads   <- .getReads(opt)
gr.reads <- gr.reads[as.logical(seqnames(gr.reads) %in% gr.contigs$qname)]
gr.contigs <- gr.contigs[gr.contigs$qname %in% unique(as.character(seqnames(gr.reads)))]

## make the trackData
suppressWarnings(cgc <- cgChain(gr.contigs, sn=gr.contigs$qname))
suppressWarnings(pac <- cgChain(gr.reads, sn=gr.reads$rname))
suppressWarnings(r2g <- gMultiply(cgc, pac, pintersect=TRUE))

suppressWarnings(tdr <- track.refgene())

## loop through and run
nams <- unique(gr.contigs$qname)
suffix = 'somatic'
sn <- as.character(seqnames(gr.reads))
if (opt$germline)
  suffix = 'germline'

grl_this = grl.sno$som
if (opt$germline)
  grl_this = grl.sno$ger

out <- mclapply(seq_along(nams), function(x) {

  print(paste("Working on contig: ", x, "of", length(nams)))
  outfile.png = file.path(opt$outdir, paste("contig_", x, '_', suffix, '.png', sep=''))
  outfile.pdf = file.path(opt$outdir, paste("contig_", x, '_', suffix, '.pdf', sep=''))

  gr.contigs.this <- gr.contigs[gr.contigs$qname == nams[x]]
  gr.reads.this <- gr.reads[sn == nams[x]]

  ## get the sequences track data
  tds <- suppressWarnings(.getTrackData(cgc, r2g, gr.contigs.this, gr.reads.this))
  tds <- c(tds, tdr)

  ## get the alignment track data
  tda <- suppressWarnings(.getTrackDataAlign(cgc, r2g, gr.contigs.this, gr.reads.this))
  
  ## set the window
  win <- streduce(gr.contigs.this)
  
  ## make the sequences plot
  png(outfile.png, width=1500, height=600)
  display(c(tda, tds), windows=win, links=grl_this, cex.ylabel=1)
  dev.off()

  pdf(outfile.pdf, width=15, height=10)
  display(c(tda, tds), windows=win, links=grl_this, cex.ylabel=1)
  dev.off()

  
  ## make the alignment plot  
  mc <- mcols(grl_this)
  mc <- mc[mc$cname == nams[x],] 
  
  ## add the figure to the plot
   expid = sprintf("BREAKPOINT_FIGURE_%s", x)
  break1 = paste(mc$chr1, ":", mc$pos1, "(", mc$strand1, "){MAPQ: ", mc$mapq1, "}", sep="")
  break2 = paste(mc$chr2, ":", mc$pos2, "(", mc$strand2, "){MAPQ: ", mc$mapq2, "}", sep="")
  lead = sprintf("Breakpoint plot of junction %s. From the top down: The red arrows at the top give the orientation of the junction, with the arrowhead pointing towards the joined DNA (away from the junction). The contig is then shown below, colored with respect to MAPQ, with white = MAPQ 0 and dark red = MAPQ 60. The contig sequence is also shown. The individual reads supporting this contig are displayed as mapped to the genome *through* their mapping to the contig. The reads were identified as potentially supporting this contig within 'snowman run', and were more thouroughly Smith-Waterman realigned in 'snowman gather'. Reads with > 6 bases on both sides of the junction are 'split' reads.", x)
  caption = paste(lead, "Breakpoint plot of junction", x, "[ Name:", mc$cname, "] -- [ Somatic breakpoint junction:", x, "] -- [ Tumor Split Reads:", mc$tsplit, "] -- [ Normal Split Reads:",
    mc$nsplit, "] -- [ Break1:", break1, "] -- [ Break2:", break2, "] -- [ Span:", mc$span,
    "] -- [ Homology:", mc$homology, "] -- [ Insertion:", mc$insertion, "] -- [ Num Times Break Found:", mc$num.dups+1, "] -- [ Number of Alignments:",
    mc$num.parts, "]")
  
  fig = newFigure(file.name(outfile.png), fileHigRes=file.name(outfile.png), exportId = expid, caption)
  return(fig);
    
}, mc.cores=opt$cores)

for (i in seq_along(out))
  r = addToResults(r, out[[i]])

# Phase 3: render report to file
print('...done with report sending to file')
outpath <- file.path(opt$outdir, paste("snowman_report", suffix, sep='_'))
writeReport(r, filename=outpath); # w/o extension

## test
if (FALSE) {
  gr.reads$tn = substring(gr.reads$rname, 1,1)
  gr.reads$rqname = gsub("[a-z]+[0-9]+_?(.*)", "\\1", gr.reads$rname)
  gr.reads$rheader = gsub("[^_]+($)", "\\1", gr.reads$rname)
  gr.reads$rqname = substring(gr.reads$rname, nchar(gr.reads$rheader)+1, nchar(gr.reads$rname))
  gr.reads$flag = as.numeric(gsub("[a-z]+([0-9]+)_?(.*)", "\\1", gr.reads$rheader))
  gr.reads.per <- split(gr.reads, seqnames(gr.reads))

  ## grab the discordant reads
  gr.test = gr.reads.per[[1]]
  tab = table(gr.test$rqname)
  gr.test.disc = gr.test[gr.test$rqname %in% names(tab[tab==2])]
  grl.test.disc = GRangesList()
  if (length(gr.test.disc))
    grl.test.disc = split(gr.test.disc, gr.test.disc$rqname)
}
