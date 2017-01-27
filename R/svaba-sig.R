#!/usr/bin/env Rscript

library(optparse)
require(VariantAnnotation)

option_list = list(
    make_option(c("-i", "--input"),        type = "character", default = "qcreport.txt",  help = "Input txt file from a snowman preprocess qcreport.txt"),
    make_option(c("-o", "--output"),       type = "character", default = "qcreport.pdf",  help = "Output pdf to generate"),
    make_option(c("-f", "--FHworkspace"),  type = "character", default = NULL,            help = "Firehose workspace to retrieve data from"),
    make_option(c("-p", "--FHpairset"),    type = "character", default = NULL,            help = "Firehose pairset to retreive data from"),
    make_option(c("-a", "--FHannotation"), type = "character", default = NULL,            help = "Firehose annoation to retreive data from")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

fhsum = sum(c(!is.null(opt$FHworkspace), !is.null(opt$FHpairset), !is.null(opt$FHannotation)))
GETFH = fhsum == 3
if (fhsum %in% c(1,2))
  print(print_help(parseobj))
  stop("If importing from FH, must specify fully the FHworkspace, FHpairset and FHannotation")
}

if (is.null(opt$input))
  stop(print_help(parseobj))

if (!file.exists(opt$input)) {
  print(print_help(parseobj))
  stop(paste("Input file does not exist", opt$input, ". Must supply path to valid qcreport.txt file (generated from snowman preprocess or snowman run"))
}

source('/home/unix/jwala/GIT/isva/Taiga/R/sourceall.R')

########################
## some functions
########################
source.all <- function() {

  suppressMessages(suppressWarnings(require(ff, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(VariantAnnotation, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(rtracklayer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(data.table, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(plyr, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(ggplot2, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(reshape2, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(GenomicRanges, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(popbio, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(bitops, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(seqinr, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Rsamtools, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(ff, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(multicore, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Biostrings, quietly=TRUE)))  
  suppressMessages(suppressWarnings(require(rtracklayer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(lattice, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(RColorBrewer, quietly=TRUE)))
  suppressMessages(suppressWarnings(require(Matrix, quietly=TRUE)))

}

round.n <- function(x, n) n * round(x / n)

lm_eqn = function(m){
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                              list(a = format(coef(m)[1], digits = 2),
                                                 b = format(coef(m)[2], digits = 2),
                                                r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

## power.law.plot
power.law.plot <- function(pspan, plotname) {
  print(paste('...making power-law plot for', plotname))
  ex <- ecdf(pspan[pspan > 0])
  #grid <- c(seq(0,10000,5), seq(10000,10^7,100))
  df = data.frame(pspan = pspan, cdf = ex(pspan), pspan.log = log(pspan,10), cdf.log = log(ex(pspan), 10))
  df <- df[df$cdf.log > -1000,] ## get rid of -Inf
  lm.s <- lm(df$cdf.log ~ df$pspan.log)
  lm.eq = lm_eqn(lm.s)

  yb = c(0.01,0.02,0.03,0.04,0.05,0.10,0.2, 0.3, 0.4,0.5,0.75,1)

  pdf(file.path(opt$outdir, plotname), width=5,height=5)
  print(g <- ggplot() + geom_point(data=df, aes(x=pspan.log, y=cdf.log), size=1, color="blue") +
    geom_smooth(data=df, aes(x=pspan.log, y=cdf.log), method="lm",se=F) +
        geom_text(aes(x=4.5,y=log(0.01,10), label=lm.eq), label=lm.eq, parse=T) +
    theme_bw() + ylab("CDF") + xlab("Span (bp)") +
    scale_y_continuous(breaks=log(yb,10), labels=format(yb,digits=2), limits=c(log(yb[1],10), 0)) +
    scale_x_continuous(breaks=seq(0,7), labels=parse(text=paste('10', seq(0,7), sep='^'))))
  dev.off()
}


##################################
################ load from data ##
##################################
if (GETFH) {
  ind <- fiss_get(opt$FHpairset, wkspace=opt$FHworkspace, type='pair')
  vcf <- ind[, opt$FHannotation]
} else {

}

print('...loading VCF files')
toload = file.exists(vcf)
cols = sample(colors(), sum(toload))
ra <- lapply(vcf[toload], function(x) {
  ab = ra_breaks(x)
  mcols(ab)$individual = basename(x)
  mcols(ab)$border = mcols(ab)$col = cols[match(x, vcf)]
  return (ab)
})

## place in data table form
ra.bound = do.call('grlbind', ra)
ra.ul <- grl.unlist(ra.bound)
ra.dt <- gr2dt(ra.ul)

## get the spans
span = ra.dt[, abs(end[1]-end[2]) * (seqnames[1] == seqnames[2]), by=grl.ix]$V1
mcols(ra.bound)$span = span
mcols(ra.bound)$span[mcols(ra.bound)$span == 0] = NA

## set the logical for methods
s  = mcols(ra.bound)$CALLER == "S"
d  = mcols(ra.bound)$CALLER == "D"
ds = mcols(ra.bound)$CALLER == "DS"

############################
## make the overlap pie chart
############################
cols = c(ovlp.hue, dran.hue, snow.hue)
levs = c("Both", "dRanger", "SnowmanSV")
df <- data.frame(Caller=c("SnowmanSV", "dRanger", "Both"), counts=c(sum(s, na.rm=T), sum(d, na.rm=T), sum(ds, na.rm=T)), stringsAsFactors=FALSE)
g <- ggplot(df, aes(x=factor(1), y=counts, fill=Caller)) + geom_bar(stat='identity', position="fill") + coord_polar(theta="y") + xlab("") + ylab("") +
  scale_fill_manual(values=cols, breaks=levs)
pdf(file.path(opt$outdir, "caller_overlap_pie.pdf")); print(g); dev.off()

############################
## make the 1/L distribution
############################
power.law.plot(span, "power_law_DS.pdf")
power.law.plot(span[which(s)], "power_law_S.pdf")
power.law.plot(span[which(d)], "power_law_D.pdf")

###############################################
## check if anything overlaps with cancer genes
###############################################
cgc.genes <- track.load('cgc')
fo <- gr.findoverlaps(ra.ul, gr.pad(cgc.genes, 5e4))
tab <- table(fo$subject.id)
gene.nums = as.numeric(names(tab))
names(tab) <- cgc.genes$gene[gene.nums]
tab.lenscale = tab / width(cgc.genes[gene.nums])

## generate PDFs for all of these
td.rg.cgc = track.refgene(genes = cgc.genes$gene, height=3)
ord = order(tab.lenscale, decreasing=T)
dir.create(file.path(opt$outdir, "CGCgenes"), showWarnings=FALSE)

dum <- sapply(seq_along(ord), function(x) 
{
  gene = names(tab.lenscale[ord[x]])
  ## find the partners
  gene.window = gr.pad(cgc.genes[cgc.genes$gene == gene], 5e4)
  rar.hits = ra.bound[ceiling(gr.findoverlaps(ra.ul, gene.window)$query.id/2)]
  windows <- streduce(grbind(gr.pad(streduce(rar.hits), 5e4), gene.window))
  nindiv = length(unique(mcols(rar.hits)$individual))
                                      
  print(paste("plotting CGC gene", gene))
  pdf(file.path(opt$outdir, "CGCgenes", paste("rank_",sprintf("%04d",x),"_",names(tab.lenscale[ord[x]]), "_NIndiv_", nindiv, ".pdf", sep="")))
  #pdf(file.path(opt$outdir, "ROS1_special.pdf"))  
  td.rg.cgc$xaxis.newline = T
  td.rg.cgc$xaxis.cex = 0.5
  td.rg.cgc$xaxis.cex.label = 0.5
  td.rg.cgc$xaxis.nticks = 2
  display(td.rg.cgc, links=rar.hits, window=windows)
  dev.off()
})

##############################
## matrix plot
##############################
grt  = gr.tile(gr.all(), w=10e6)
mat <- sig.tri.matrix(ra.bound, grt, log=TRUE)
td.mat  <- do.call('trackData', c(allopts, list(grt, mdata=mat, triangle=TRUE,
                  cmap.min=min(mat), cmap.max=max(mat)+0.1, track.name='Triangle',
                  height=25, sep.lwd=0.5, m.bg.col='white',
                  track.name='Breakpoint-pair Heatmap', islog=TRUE, xaxis.nticks=0,
                  xaxis.prefix="", xaxis.chronly=TRUE)))
pdf(file.path(opt$outdir, "overview_matrix.pdf"), width=12, height=12)
display(td.mat, windows=streduce(gr.all()))
dev.off()

## per chrom
grt  = gr.tile(gr.all(), w=10e6)
mat <- sig.tri.matrix(ra.bound, grt, log=TRUE)
td.mat  <- do.call('trackData', c(allopts, list(grt, mdata=mat, triangle=TRUE,
                  cmap.min=min(mat)+0.1, cmap.max=max(mat)+0.1, track.name='Triangle',
                  height=25, sep.lwd=0.5, m.bg.col='white',
                  track.name='Breakpoint-pair Heatmap', islog=TRUE, xaxis.nticks=0,
                  xaxis.prefix="", xaxis.chronly=TRUE)))
dum <- lapply(seq(23), function(x) {
  print(paste("plotting for chr", x))
  pdf(file.path(opt$outdir, paste("chr",x,"_matrix.pdf",sep="")), width=12, height=12)
  display(td.mat, windows=streduce(gr.all()[x]))
  dev.off()
})


################################################
## Start the significance testing ##
################################################
## load data for 1e6 binned mapping
dload <- sig.load.bins('1e6')
grt <- dload$grt
mat.background = sig.2d.map.background(dt2gr(grt), nomap=FALSE)

fo <- sig.bin.counts(ra.bound, dt2gr(grt)) ## place events into bins
fo[, sample.counts := nrow(.SD), by='sample'] ## get total number of binned-events per sample
setkey(fo, sample) ## sort by sample count
M <- structure(as.numeric(fo$sample.counts[!duplicated(fo)]), names=fo$sample[!duplicated(fo)]) ## M is number mutations per sample
M <- M[order(names(M))]
exclude <- names(M[M < 1]) ## drop really weak samples
M <- M[!names(M) %in% exclude]
ra.bound <- ra.bound[!mcols(ra.bound)$individual %in% exclude]

## normalize background probablity matrix by total intra/inter chrom rate
empirical.background.prob <- sig.generate.normalized.matrices(ra.bound, dload$mat, M, single.mat=TRUE)
## get the total vector of background F-scores
system.time(f.back <- sig.2d.permute(empirical.background.prob, M, num.permute = 1000, method='1d'))
f.back2 <- f.back[f.back > 0]
## plot the distribution of f scores
f.back2 <- f.back[f.back > 0]
pdf(file.path(opt$outdir, "fscore_distribution.pdf"), height=3.5, width=10)
hist(f.back2, col='black', breaks=300, xlab='F-score', ylab='Count', main='F-score distribution under NULL')
dev.off()

## get a trackData for the bin probabilities
#all.bin.probs <- colSums(do.call('rbind', lapply(all.mats, colSums))) / length(all.mats)
all.bin.probs <- colSums(empirical.background.prob)
grt$binprobs <- all.bin.probs
grt$col <- grt$border <- chr_colors[as.character(seqnames(dt2gr(grt)))]
grt$binprobs <- -log(grt$binprobs, 10)
td.binprob <- trackData(dt2gr(grt), y.field='binprobs', xaxis.newline=FALSE, xaxis.chronly=TRUE, track.name='-Log(Bin Background Probability)',
                        y0=3.3, y1=3.9, circles=TRUE, lwd.border=0.5, xaxis.nticks = 0, xaxis.prefix="")
pdf(file.path(opt$outdir, "td_bin_background_probabiltiies.pdf"), width=10)
display(td.binprob, windows=gr.all())
dev.off()

## get the actual f-score
f.real <- sig.fscore(grl, dt2gr(grt), empirical.background.prob)

## get the p and q values
pq = sig.fscore2pval(f.back, f.real)
grt$pval <- pq$pval
grt$qval <- pq$qval

## make the qq-plot
pdf(file.path(opt$outdir, "qqplot.pdf"), width=5, height=5)
qq_pval(grt$pval)
dev.off()

grt$n.log10.q <- -log(grt$qval, 10)
grt$col <- grt$border <- chr_colors[as.character(seqnames(dt2gr(grt)))]
td.sig <- trackData(dt2gr(grt), y.field='n.log10.q', circles=TRUE, track.name='-log10(Q-value)', lwd.border=0.5, xaxis.chronly=T, y0=0, xaxis.prefix="", xaxis.nticks=0)
ppdf(display(td.sig, windows=gr.all()))

td.refgenes <- track.refgene()
dir.create(file.path(opt$outdir, "topsighits"), showWarnings=FALSE)
### plot all the hits above -log10(q) = 1
which.hits <- which(grt$n.log10.q) > 1
dum <- lapply(seq_along(which.hits), function(x)
{

  y = which.hits[x]
  window = dt2gr(grt)[y]
  genes <- unique(gr.genes$gene[gr.findoverlaps(gr.genes, window)$query.id])
  genes.string <- paste(genes, collapse = " ")
  print(paste("Plotting top significance bin", x, "of", sum(grt$n.log10.q > 1), "which is region", grt$seqnames[y],":", grt$start[y], "-", grt$end[y], "which contains genes", genes.string))
  
  rar.hits <- ra.bound[ceiling(gr.findoverlaps(ra.ul, window)$query.id/2)]
  disp.window <- streduce(grbind(gr.pad(streduce(rar.hits),5e4),window))

  pdf(file.path(opt$outdir, "topsighits", paste("hit_", x, "q_",sprintf("%04f",grt$n.log10.q[y]),".pdf", sep="")))
  td.rg.cgc$xaxis.newline = T
  td.rg.cgc$xaxis.cex = 0.5
  td.rg.cgc$xaxis.cex.label = 0.5
  td.rg.cgc$xaxis.nticks = 2
  td.rg.cgc$cex.tick = 0.5
  td.rg.cgc$track.name = "Cancer Genes"
  td.refgenes$track.name = "All Genes"
  display(c(td.rg.cgc, td.refgenes), links=rar.hits, window=disp.window)
  dev.off()
  hitnum = hitnum + 1
})
