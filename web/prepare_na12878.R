## function for analyzing na12878 deletions
.na12878_olap <- function(x) {
  
  # if ("indel" %in% names(x)) {
  #   
  #   x2 <- x$indel[x$indel$SPAN > 50 & x$indel$ETYPE == "del"]
  #   x3_1 <- GRanges(seqnames(x2), IRanges(start(x2), width=1), strand="+", id=seq_along(x2))
  #   x3_2 <- GRanges(seqnames(x2), IRanges(end(x2), width=1), strand="-", id=seq_along(x2))
  #   b <- c(x3_1, x3_2)
  #   b <- split(b, b$id)
  #   mcols(b)$deltype <- TRUE
  #   mcols(b)$SPAN <- x2$SPAN
  #   mcols(b)$pacbio <- FALSE
  #   x$grl <- grlbind(x$grl, b)
  #   x$dt$id = x$dt$SCTG
  #     
  # }
  
  if ("deltype" %in% colnames(mcols(x$grl))) {
    dd = x$grl[mcols(x$grl)$deltype]
  } else if ("SVTYPE" %in% colnames(x$grl)) { ## delly
    dd = x$grl[mcols(x$grl)$SVTYPE=="DEL"]
  } else {
    dd = x$grl
    x$dt$deltype <- TRUE
    mcols(dd)$deltype <- TRUE
  }
  
  if ("SVLEN" %in% colnames(mcols(dd))) {
    mcols(dd)$SPAN = mcols(dd)$SVLEN
    x$dt[x$dt$ALT >= 6]
    dd <- dd[!is.na(mcols(dd)$ALT) & mcols(dd)$ALT >= 6]
  }
  
  if (!"id" %in% colnames(x$dt))
    x$dt$id <- seq(nrow(x$dt))
  
  if ("pacbio" %in% colnames(mcols(dd)))
    pacbio <- mcols(dd)$pacbio
  else
    pacbio <- rep(FALSE, length(dd))
  print(paste("Pacbio hits", sum(pacbio)))
  
  suppressWarnings(ro2 <- ra.overlaps(dd, truth.NA12878$grl, pad=2e2, ignore.strand=TRUE))
  TP2e <- dd[unique(c(ro2[,'ra1.ix'], which(pacbio)))]
  TP2 <- length(unique(ro2[,'ra2.ix'])) ## TRUE POSITIVE
  FP2e <- dd[setdiff(seq_along(dd), c(ro2[,'ra1.ix'], which(pacbio)))]
  FP2 = sum(mcols(FP2e)$SPAN >= 100 & mcols(FP2e)$SPAN < 1e5)
  
  TP1=FP1=FP1e=TP1e=0
  
  #setkey(x$dt, seqnames, start)
  return(list(dels=x$dt[!duplicated(x$dt$id) & x$dt$deltype & x$dt$SPAN >= 50], dels.grl=x$grl[mcols(x$grl)$deltype], FP1=FP1, FP2=FP2, TP1=TP1, TP2=TP2, ro1=NULL, ro2=ro2, TP1e=TP1e,TP2e=TP2e, FP1e=FP1e,FP2e=FP2e))
}

#################################
### TRUTH SET FROM LUMPY PAPER
#################################



## Load the truth set from the LUMPY paper, the set with the added calls validated by PacBio/Moleculo
ff <- fread("data/lumpy_na12878_truthset4.bedpe")
setnames(ff, paste0("V",seq(10)), c("seqnames", "start1","end1", 'altchr','altpos','altend',"ID","SPAN","strand","altstrand")) 
ff[ , seqnames := gsub("chr","", seqnames)]
ff[ , altchr := gsub("chr","", altchr)]
ff2 <- data.table::copy(ff) 
setnames(ff2,c("seqnames","start1","end1","strand","altchr","altpos","altend","altstrand"),c("altchr","altpos","altend", "altstrand","seqnames", "start1","end1","strand"))
ff2$strand <- "-"
ff_double <- rbind(ff,ff2)
ff_double[, start := (start1+end1)/2, by=ID]
ff_double[, end := start]

dels_double <- sort(gr.fix(gr.nochr(dt2gr(ff_double)), si))
grl.ff <- split(dels_double, dels_double$ID)

saveRDS(list(grl=grl.ff, dt=ff), "data/lumpy_na12878_truthset4.rds")
truth.NA12878 <- readRDS("data/lumpy_na12878_truthset4.rds")


#################
#### LUMPY
#################

## load the LUMPY calls as generated from FH
lumpy <- load_lumpy("data/lumpy_f24_na12878.vcf")[[1]]
## limit to +- orientation (del type)
ll <- gr2dt(grl.unlist(lumpy))
ll[, c("id.1","PE.1","IMPRECISE.1","SECONDARY.1","SR.1") := NULL]
ll[, deltype := strand[1] == '+' && strand[2] == "-" && seqnames[1] == seqnames[2], by=grl.ix]
setkey(ll, grl.ix)
gr <- dt2gr(ll)
mcols(gr) <- NULL
grl <- split(gr, ll$grl.ix)
mcols(grl) <- ll[!duplicated(grl.ix), .(id, PE, IMPRECISE, SECONDARY, SR, SPAN, deltype)]
lum<-list(grl=grl, dt=ll)

## make the LUMPY ROC
lum.out <- rbindlist(lapply(seq(2, 20, 1), function(x) {
  print(x)
  
  llo <- lum
  llo$grl <- lum$grl[ix <- mcols(lum$grl)$SR + mcols(lum$grl)$PE >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="LUMPY"))
  return(dt)
  
}))
saveRDS(lum.out, "data/lumpy_roc_na12878.rds")

#################
## PINDEL
#################
ff <- fread("data/pindel.na12878.raw")
ff[, SUM_MS := as.integer(gsub("SUM_MS ([0-9]+)", "\\1", V16))]                                                                                                                                                                                                               
ff[, ALT := as.integer(gsub("Supports ([0-9]+)", "\\1", V9))]                                                                                                                                                                                                                 
ff[, MEAN_MAPQ := ifelse(ALT > 0, SUM_MS / ALT, 0)]                                                                                                                                                                                                                           
ff[, SVTYPE := gsub("([A-Z]+).*","\\1",V2)]                                                                                                                                                                                                                                   
ff[, seqnames := gsub("ChrID (.*?)","\\1",V4)]                                                                                                                                                                                                                                
ff[, SPAN := as.integer(gsub(".*?([0-9]+)","\\1",V2))]                                                                                                                                                                                                                        
ff[, start := as.integer(gsub("BP ([0-9]+)", "\\1", V5))]                                                                                                                                                                                                                     
ff[, NT_LEN := as.numeric(gsub("NT ([0-9]+).*", "\\1", V3))]                   

pindel.dels <- ff[SPAN > 30 & SVTYPE=="D" & MEAN_MAPQ >= 20 & ALT >= 2 & NT_LEN < 10]
gr <- with(pindel.dels, GRanges(c(seqnames, seqnames), IRanges(c(start, V6), width=1), strand=rep(c("+","-"), each=nrow(pindel.dels)), id=rep(seq(nrow(pindel.dels)),2))) 
grl.pindel <- split(gr, gr$id) 
mcols(grl.pindel) <- pindel.dels[,.(MEAN_MAPQ, ALT, SVTYPE, SPAN, NT_LEN)]
ss=list(dt=pindel.dels, grl=grl.pindel)

ss.out <- rbindlist(lapply(2:20, function(x) {
  print(x)
  
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$ALT) & mcols(ss$grl)$ALT >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Pindel"))
}))
saveRDS(ss.out, "data/pindel_roc_na12878.rds")

###################
### DELLLY
###################
delly <- load_delly("data/delly.na12878.vcf")[[1]]
ss = list(grl=delly, gr2dt(grl.unlist(delly)))

ss$grl <- ss$grl[mcols(ss$grl)$SVTYPE =="DEL" & mcols(ss$grl)$SPAN > 50]
ss.out <- rbindlist(lapply(2:20, function(x) {
  print(x)
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$DISC) & (mcols(ss$grl)$DISC + mcols(ss$grl)$SPLIT) >= x]
  ll2 <- .na12878_olap(llo)
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="DELLY (RP + SR)"))
  
  llo <- ss
  llo$grl <- ss$grl[ix <- !is.na(mcols(ss$grl)$DISC) & mcols(ss$grl)$DISC >= x]
  ll2 <- .na12878_olap(llo)
  dt <- rbind(dt, with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="DELLY (RP)")))
  
}))
saveRDS(ss.out, "data/delly_roc_na12878.rds")


###################
### SNOWMAN
###################
## SNOWMAN

snow <- load_snowman("data/snowman.na12878.sv.vcf")[[1]]
snowi <- load_indel("data/snowman.na12878.indel.vcf")[[1]]

## get just the dels                                                                                                                                                                                                                                                          
rr <- gr2dt(grl.unlist(snow)) 
rr[, deltype := strand[1] == '+' && strand[2]=='-' && seqnames[1] == seqnames[2] && nchar(INSERTION) < SPAN, by=grl.ix]

mcols(snow)$deltype <- FALSE                                                                                                                                                                                                                                                  
mcols(snow)$deltype[rr$deltype[!duplicated(rr$grl.ix)]] <- TRUE                                                                                                                                                                                                               
snow <- snow[mcols(snow)$NORMALT >= 2 & mcols(snow)$NORMCOV > 10] # & mcols(snow)$NORMALT / mcols(snow)$NORMCOV > 0.2]                                                                                                                                                        
#g <- flag.plot(xsv=snow[rr$deltype[!duplicated(rr$grl.ix)]], xindel=snowi, e=dels, pad=2e3)                                                                                                                                                                                  

s = gr2dt(grl.unlist(snow))                                                                                                                                                                                                                                                   
s[, altchr := rev(seqnames), by=grl.ix]                                                                                                                                                                                                                                       
s[, altpos := rev(start), by=grl.ix]                                                                                                                                                                                                                                          
s[, altstrand := rev(strand), by=grl.ix]                                                                                                                                                                                                                                      
s2 <- gr2dt(snowi)                                                                                                                                                                                                                                                            
s2[, strand := '+']                                                                                                                                                                                                                                                           
s2[, altstrand := '-']                                                                                                                                                                                                                                                        
s2[, altchr := seqnames]                                                                                                                                                                                                                                                      
s2[, altpos := end]                                                                                                                                                                                                                                                           
s2[, deltype := ETYPE == 'del']                                                                                                                                                                                                                                               
b=intersect(colnames(s2), colnames(s))                                                                                                                                                                                                                                        
ss <- rbind(s2[,b, with=FALSE],s[,b, with=FALSE])                                                                                                                                                                                                                             
print("saving")                                                                                                                                                                                                                                                               

x2 <- snowi[snowi$SPAN > 50 & snowi$ETYPE == "del"]                                                                                                                                                                                                                           
x3_1 <- GRanges(seqnames(x2), IRanges(start(x2), width=1), strand="+", id=seq_along(x2))                                                                                                                                                                                      
x3_2 <- GRanges(seqnames(x2), IRanges(end(x2), width=1), strand="-", id=seq_along(x2))                                                                                                                                                                                        
b <- c(x3_1, x3_2)                                                                                                                                                                                                                                                            
b <- split(b, b$id)                                                                                                                                                                                                                                                           
mcols(b)$deltype <- TRUE                                                                                                                                                                                                                                                      
mcols(b)$SPAN <- x2$SPAN                                                                                                                                                                                                                                                      
grl.snowi <- b                                                                                                                                                                                                                                                                
mcols(grl.snowi)$EVDNC <- "INDEL"                                                                                                                                                                                                                                             
mcols(grl.snowi)$NORMALT <- x2$AD                                                                                                                                                                                                                                             
saveRDS(ss <- list(dt=ss[,.(seqnames, start, strand, altchr, altpos, altstrand, SPAN, deltype, NM, SCTG, MAPQ, REPSEQ)], grl=grlbind(snow, grl.snowi), indel=snowi), "data/snowman_na12878.rds") 

###############
## SNOWMAN
###############
ss <- readRDS("data/snowman_na12878.rds")
##ss$grl <- ss$grl[-unique(c(which(pmax(mcols(ss$grl)$DISC_MAPQ, mcols(ss$grl)$MAPQ) <= 10),which(pmax(mcols(ss$grl)$NM, mcols(ss$grl)$MATENM) >= 10)))]
ss.out <- rbindlist(lapply(seq(2,20,1), function(x) {
  print(x)
  
  llo <- ss
  llo$grl <- ss$grl[ix <- mcols(ss$grl)$NORMALT >= x & mcols(ss$grl)$EVDNC %in% c("ASSMB", "ASDIS")]
  #llo$indel <- ss$indel[ss$indel$DP >= x]
  ll2 <- .na12878_olap(llo)
  
  llo <- ss
  llo$grl <- ss$grl[ix <- mcols(ss$grl)$NORMALT >= x]
  #llo$indel <- ss$indel[ss$indel$DP >= x]
  ll3 <- .na12878_olap(llo)
  
  dt <- with(ll2, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Snowman (AS)"))
  dt <- rbind(dt, with(ll3, data.table(x=x, TP1=TP1, TP2=TP2, FP1=FP1,FP2=FP2, caller="Snowman (RP + AS)")))
  return(dt)
}))
saveRDS(ss.out, "data/snowman_roc_na12878.rds")
