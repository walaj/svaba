source("batch.functions.R")
library(skitools)
library(gUtils)
library(gTrack)
library(data.table)
library(skidb)

Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")

## load the NA12878 data
truth.NA12878       <- readRDS("data/lumpy_na12878_truthset4.rds")
snowman_roc_na12878 <- readRDS("data/snowman_roc_na12878.rds")
pindel_roc_na12878  <- readRDS("data/pindel_roc_na12878.rds")
delly_roc_na12878   <- readRDS("data/delly_roc_na12878.rds")
lumpy_roc_na12878   <- readRDS("data/lumpy_roc_na12878.rds")

#ggplot(data=rbindlist(list(delly_roc_na12878, lumpy_roc_na12878, pindel_roc_na12878, snowman_roc_na12878)), aes(x=FP2, y=TP2, color=caller)) + geom_point() + geom_line() + theme_bw() + xlab("Not in validation set") + scale_x_continuous(breaks=seq(0, 3000, by=250)) + ylab("True Positive") + scale_y_continuous(breaks=seq(1200,2900,by=100)) + coord_cartesian(xlim=c(0,3000), ylim=c(1500,2900))

###########################
## SIMLULATED
###########################

print("...loading simulated data")
events.d1 <- fread("data/events.d1.txt")
setnames(events.d1, paste0("V",seq(11)), c("chr","pos","altchr","altpos","strand","altstrand","dummy","span","class","ins_seq","ID"))
events.d1[class == "del"]$ins_seq <- ""
events.d1[class == "ins"]$ins_seq <- events.d1[class == "ins"]$dummy
events.d1[, dummy := NULL]
setkey(events.d1, chr, pos)
gr.events <- sort(gr.fix(with(events.d1, GRanges(c(chr,altchr), IRanges(c(pos,altpos), width=1), strand=c(strand,altstrand), id=c(ID,ID), span=c(span,span), type=c(class,class), ins_seq=c(ins_seq, ins_seq))),si))

#.prepare_bps(2)
#.prepare_bps(5)
#.prepare_bps(10)

## load pindel
#pindel <- readRDS("data/pindel_sim.rds")

snow_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr=dat
  datr$sv = datr$sv[mcols(datr$sv)$TUMALT >= x]
  datr$indel = datr$indel[mcols(datr$indel)$AD >= x]
  dt <- flag.plot(xsv=datr$sv, xindel=datr$indel, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Snowman"
  dt$ALT = x
  return(dt)
})
rb <- rbindlist(snow_sim_results)

ff <- readRDS("data/pindel_sim.rds")[[2]]
pindel_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[ff$TUMALT >= x & ff$NORMAL <= 0]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Pindel"
  dt$ALT = x
  return(dt)
})
pindel_sim_results <- rbindlist(pindel_sim_results)

ff <- readRDS("data/platypus_sim.rds")[[2]]
platypus_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[ff$NV_t >= x]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Platypus"
  dt$ALT = x
  return(dt)
})
platypus_sim_results <- rbindlist(platypus_sim_results)

ff <- readRDS("data/lumpy_sim.rds")[[2]]
mcols(ff)$TUMALT <- mcols(ff)$SR + mcols(ff)$PE
mcols(ff)$SPAN <- mcols(ff)$span
lumpy_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[mcols(ff)$TUMALT >= x]
  dt <- flag.plot(xsv=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Lumpy"
  dt$ALT = x
  return(dt)
})
lumpy_sim_results <- rbindlist(lumpy_sim_results)

ff <- readRDS("data/strelka_sim.rds")[[2]]
ff$SPAN <- ff$CALC_SPAN
strelka_sim_results <- lapply(seq(2,10,1), function(x) {
  print(x)
  datr= ff[mcols(ff)$TALT >= x]
  dt <- flag.plot(xindel=datr, e=gr.events, pad=10, fname="noplot")
  dt$caller = "Strelka"
  dt$ALT = x
  return(dt)
})
strelka_sim_results <- rbindlist(strelka_sim_results)

ggplot(data=rbindlist(list(rb[Group=="SV"],lumpy_sim_results[Group=="SV"], pindel_sim_results[Group=="SV"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + scale_x_continuous(breaks=seq(0, 15, by=1)) + coord_cartesian(xlim=c(0,15))
ggplot(data=rbindlist(list(rb[Group=="Indel"],strelka_sim_results[Group=="Indel"],pindel_sim_results[Group=="Indel" & FP < 1000],platypus_sim_results[Group=="Indel"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + ylab("True Positive") 
ggplot(data=rbindlist(list(rb[Group=="Medium"],lumpy_sim_results[Group=="Medium"], pindel_sim_results[Group=="Medium" & FP < 500],platypus_sim_results[Group=="Medium"])), aes(x=FP, y=TP, color=caller))+ geom_point() + geom_line() + theme_bw() + xlab("False Positive") + ylab("True Positive") # + scale_y_continuous(breaks=seq(0,3000,500)) + scale_x_continuous(breaks=seq(0,5,20))

print("done loading data")

