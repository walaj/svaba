#!/usr/bin/env Rscript

library(optparse)

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = "qcreport.txt",  help = "Input txt file from a snowman preprocess qcreport.txt"),
    make_option(c("-o", "--output"), type = "character", default = "qcreport.pdf",  help = "Output pdf to generate")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

if (!file.exists(opt$input)) {
  print(print_help(parseobj))
  stop(paste("Input file does not exist", opt$input, ". Must supply path to valid qcreport.txt file (generated from snowman preprocess or snowman run"))
}

print(opt)

require(ggplot2)
require(reshape2)
require(gridExtra)

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

colnames(df.mapq) <- seq(from=0,to=60)
colnames(df.nm) <- seq(from=0,to=ncol(df.nm)-1)
colnames(df.isize) <- seq(from=0,to=ncol(df.isize)-1)
colnames(df.xp) <- seq(from=0, to=ncol(df.xp)-1)
colnames(df.as) <- seq(from=0, to=ncol(df.as)-1)
colnames(df.len) <- seq(from=0, to=ncol(df.len)-1)
colnames(df.phred) <- seq(from=0, to=ncol(df.phred)-1)
colnames(df.clip) <- seq(from=0, to=ncol(df.clip)-1)
readg <- sapply(rg, function(x) x$readgroup)

df.mapq$readgroup  <- readg
df.nm$readgroup    <- readg
df.isize$readgroup <- readg
df.xp$readgroup    <- readg
df.as$readgroup    <- readg
df.phred$readgroup <- readg
df.len$readgroup   <- readg
df.clip$readgroup  <- readg

g.mapq  <- ggplot(df.mapq.m  <- melt(df.mapq, id='readgroup'),  aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(30,61))    + ylab('Reads') + xlab('Mapping Quality') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.nm    <- ggplot(df.nm.m    <- melt(df.nm, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,50))     + ylab('Reads') + xlab('NM Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.isize <- ggplot(df.isize.m <- melt(df.isize, id='readgroup'), aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(10,2001))  + ylab('Reads') + xlab('InsertSize') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.xp    <- ggplot(df.xp.m    <- melt(df.xp, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,100))    + ylab('Reads') + xlab('XP Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.as    <- ggplot(df.as.m    <- melt(df.as, id='readgroup'),    aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,100))    + ylab('Reads') + xlab('AS Tag') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.len   <- ggplot(df.len.m   <- melt(df.len, id='readgroup'),   aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(20,102))   + ylab('Reads') + xlab('Read Length') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.clip  <- ggplot(df.clip.m  <- melt(df.clip, id='readgroup'),  aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,85))     + ylab('Reads') + xlab('Clipped bases') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))
g.phred <- ggplot(df.phred.m <- melt(df.phred, id='readgroup'), aes(x=as.numeric(variable), y=pmax(log(value,10),0), group=readgroup, color=readgroup)) + geom_line() + theme_bw() + scale_x_continuous(limits=c(0,43))     + ylab('Reads') + xlab('Mean read Phred quality') + scale_y_continuous(breaks=seq(0,8), labels=parse(text=paste('10', seq(0,8), sep='^'))) + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))

df.isize.m$variable <- as.numeric(as.character(df.isize.m$variable))
df.disco <- df.isize.m[df.isize.m$variable > 0 & df.isize.m$variable < 800, ]
g.isize <- ggplot(df.disco, aes(x=variable, y=value, color=readgroup)) + geom_line() + theme_bw() + ylab('Reads') + xlab('Insert Size') + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9))

## get percentages
# disc, proper
df.disc.pie <- data.frame(rbind(rowSums(df.isize[, (800+1):(2000+1)]), rowSums(df.isize[, 1:800])), class=c("\"Discordant\" (I > 800)", "\"Proper\"(I < 800)"))
colnames(df.disc.pie) <- c(as.character(levels(readg)), "Class")
df.disc.pie.m <- melt(df.disc.pie, id='Class')
g.disc.pie <- ggplot(df.disc.pie.m, aes(x=factor(1), y=value, fill=Class)) + geom_bar(stat='identity', position="fill") + facet_wrap( ~ variable) + coord_polar(theta="y") + xlab("") + ylab("")  + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9), legend.position="bottom")

df.clip.pie <- data.frame(rbind(rowSums(df.clip[, 6:102]), rowSums(df.clip[, 1:5])), class=c("\"Clipped\" (clip >= 5)", "\"Matched\"(clip < 5)"))
colnames(df.clip.pie) <- c(as.character(levels(readg)), "Class")
df.clip.pie.m <- melt(df.clip.pie, id='Class')
g.clip.pie <- ggplot(df.clip.pie.m, aes(x=factor(1), y=value, fill=Class)) + geom_bar(stat='identity', position="fill") + facet_wrap( ~ variable) + coord_polar(theta="y") + xlab("") + ylab("")  + theme(legend.text = element_text(size=6), legend.title = element_text(size=6), text = element_text(size=9), legend.position="bottom")

pdf(opt$output, width=22, height=12)
print(grid.arrange(g.mapq, g.nm, g.isize, g.len, g.clip, g.phred, g.disc.pie, g.clip.pie, ncol=3))
dev.off()

