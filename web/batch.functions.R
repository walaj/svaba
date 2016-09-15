load_discovar <- function() {

  ff <- fread("/xchip/gistic/Jeremiah/Projects/Discovar/Assembly2VCF/Feb09/assembly.bps.txt")

  ff[, SL_t := as.numeric(gsub(".*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:(.*?)", "\\1", V35))]
  ff[, SL_n := as.numeric(gsub(".*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:(.*?)", "\\1", V34))]
  ff[, TUMREF := as.numeric(gsub("0/1:.*?:(.*?):.*", "\\1", V35))]
  ff[, NORMREF := as.numeric(gsub("0/1:.*?:(.*?):.*", "\\1", V34))]
  ff[, TUMALT := as.numeric(gsub("0/1:(.*?):.*", "\\1", V35))]
  ff[, NORMALT := as.numeric(gsub("0/1:(.*?):.*", "\\1", V34))]

  ## remove exact duplicate breakpoints
  setkey(ff, V1, V2, V4, V5)
  ff <- unique(ff)
  
  ## make the span
  ff[, span := ifelse(V1==V4,abs(V5-V2), 1e9)]

  ## convert to GRanges
  gr.ff <- with(ff, GRanges(V1, IRanges(V2, width=1),SPAN=span, EVDNC=V23, CONF=V22, SOM=ff$V26, REP=V30, SCTG=V20))

  ##
  #gr.ff_som <- with(ff[V26 == 1 & (V22 == "PASS" | V22 == "NOLOCAL")], GRanges(V1, IRanges(V2, width=1), REP=V30, SPAN=span, SCTG=V20, EVDNC=V23))
  som_contigs <- ff[V26 == 1 & TUMALT >= 3 & NORMALT <= 1 & (V22 == "PASS" | V22 == "NOLOCAL"), unique(V20)] ## get contig IDs
  setkey(ff, V20)
  gr.ff_som <- sort(gr.fix(with(ff[J(som_contigs)], GRanges(c(V1, V4), IRanges(c(V2,V5), width=1), strand=c(V3,V6), SCTG=c(V20,V20), SPAN=c(span,span), EVDNC=c(V23, V23), TUMALT=c(TUMALT, TUMALT), NORMALT=c(NORMALT, NORMALT), TUMREF=c(TUMREF,TUMREF), NORMREF=c(NORMREF,NORMREF))), si))
  gr.ff_som <- gr.ff_som[gr.ff_som$NORMALT <= 1]
  #gr.ff_som <- with(ff[ff$V26 == 1 & (ff$V22 == "PASS")], GRanges(V1, IRanges(V2, width=1), REP=V30, SPAN=span, SCTG=V20, EVDNC=V23))  

  ## remove somatic events that overlap with germline
  fo <- gr.findoverlaps(gr.ff[gr.ff$CONF=="PASS" & gr.ff$SOM==0 & gr.ff$EVDNC=="INDEL"], gr.ff_som+5)
  keep <- gr.ff_som[setdiff(seq_along(gr.ff_som), fo$subject.id)]
  filtered_out <- gr.ff_som[gr.ff_som$SCTG %in% gr.ff_som$SCTG[fo$subject.id]]
  gr.ff_som <- keep

  return (list(all=gr.ff, som=gr.ff_som))
}

# overlap a series of calls
venn <- function(events, pad=10) {

  if (is.null(names(events)))
    names(events) <- LETTERS[seq_along(events)]
  
  gr <- do.call('grbind', events)
  gr$method <- rep(names(events), sapply(events, length))
  
  agreements <- ovl <- list()
  
  for (i in names(events))
    {
      print(paste('...working on', i))
      noti <- gr[gr$method != i]
      grt <- gr[gr$method == i]

      fo <- gr.findoverlaps(grt, noti + pad, ignore.strand=TRUE)

      fo$method <- noti$method[fo$subject.id]
      ## make sure one event doesn't overlap with two from same method
      fo <- fo[!duplicated(paste(fo$method, fo$query.id))]

      ## convert to data.table for ops
      fo <- data.table(seqnames=as.character(seqnames(fo)), start=start(fo), end=end(fo), query.id=fo$query.id, method=fo$method)
      ## get every overlap
      fo[, methods := paste(sort(method), collapse="_"), by=query.id]

      ## number with 1, 2, 3, etc agreements
      ac <- c("0"=length(setdiff(seq_along(grt), fo$query.id)), table(table(fo$query.id)), table(fo$method), table(fo$methods[!duplicated(fo$query.id)]))
      suppressWarnings(agreements[[i]] <- ac)

      ## trim it down
      fo[, method := NULL]
      setkey(fo, query.id)

      grt$methods <- ""
      grt$methods[fo$query.id] <- fo$methods
      suppressWarnings(ovl[[i]] <- grt)

      
      
    }
  
  return(list(a=agreements,fo=ovl))

}



gr2dt <- function(gr, basic=FALSE) {
  if (any(class(gr)=='data.table'))
    return(gr)
  out <- with(gr, data.table(seqnames=as.character(seqnames(gr)),
                            start=start(gr), end=end(gr), strand=as.character(strand(gr))))
  if (!basic && ncol(mcols(gr)))
    out <- cbind(out, as.data.frame(mcols(gr)))
  return(out)
}


## load the lumpy outputs
load_lumpy <- function(files) {

  lumpy <- lapply(files, function(x) {
    
    print(basename(x))
    if (!file.exists(x)) {
      print(paste("File does not exist",x))
      return (GRangesList())
    }
    dvcf <-   readVcf(x, "hg19")
    if (exists("rowRanges"))
      d <- rowRanges(dvcf)
    else
      d <- rowData(dvcf)    
    mcols(d) = cbind(mcols(d), info(dvcf))

    if (ncol(geno(dvcf)[[2]])==2) {
      somatic <- geno(dvcf)[[2]][,2] == 0
      d <- d[somatic]
    }
    if (length(d) == 0)
      return (GRangesList())
    ALT = unlist(d$ALT)
    bnd = d$SVTYPE=="BND"
    lbe <- GRanges(seqnames(d[bnd]), IRanges(start(d[bnd]), width=1), id = gsub("([0-9]+)_[0-9]", '\\1', unlist(d$MATEID)), PE=unlist(d$PE)[bnd], IMPRECISE=d$IMPRECISE[bnd], SECONDARY=d$SECONDARY[bnd], SR=unlist(d$SR[bnd]), strand=ifelse(grepl("^\\[", ALT[bnd]) | grepl("^\\]", ALT[bnd]), '-', '+'))
    lnb <- GRanges(rep(seqnames(d[!bnd]), 2), IRanges(c(start(d[!bnd]), mcols(d[!bnd])$END),width=1), id=paste0("A",rep(seq(sum(!bnd)), 2)), PE=unlist(d$PE)[!bnd], IMPRECISE=d$IMPRECISE[!bnd], SECONDARY=d$SECONDARY[!bnd], SR=unlist(d$SR[!bnd]), strand=c(ifelse(ALT[!bnd]=="<DUP>", '-','+'),ifelse(ALT[!bnd]=="<DUP>", '+','-')))
    
    l <- c(lbe, lnb)
    dt <- data.table(id=l$id, start=start(l), end=end(l), chr=as.character(seqnames(l)), PE=l$PE, IMPRECISE=l$IMPRECISE, SECONDARY=l$SECONDARY, SR=l$SR, strand=as.character(strand(l)))
    dt[, SPAN := ifelse(chr[1] != chr[2], as.integer(1e9), abs(start[1] - end[2])), by=id]
    setkey(dt, id)
    dt[, altchr := rev(chr), by=id]
    dt[, altpos := rev(start), by=id]
    dt[, altstarnd := rev(strand), by=id]
    dt <- unique(dt)
    
    l <- split(l, l$id)
    mcols(l) <- dt[,.(id, PE, IMPRECISE, SECONDARY, SR, SPAN)]
#    mcols(l)$id <- unique(grl.unlist(l)$id)
#    mcols(l)$SPAN = dt[mcols(l)$id]$span
    return (l)
    
  })
  
  return (lumpy)
  
}

load_snowman <- function(files, mc.cores=mc.cores, unlist=FALSE, bad.remove=TRUE) {

  snow4 <- lapply(files, function(x) {
    if (!file.exists(x)) {
      print(paste("File does not exist",x))
      if (unlist)
        return (GRanges())
      else
        return (GRangesList())
    }
    
    #print(paste(basename(x), "is", match(x, files), "of", length(files)))
    #aa <- skitools::ra_breaks(x)
    #aa2 <- gr.flipstrand(grl.unlist(aa))
    #dd <- mcols(aa)
    #grlix <- aa2$grl.ix
    #mcols(aa2) <- NULL
    #aa <- split(aa2, grlix)
    #mcols(aa) <- dd

    ## dont continue if empty
    if (!as.numeric(system(paste("grep -v '^#'", x, "| wc -l"), intern=TRUE))) {
      if (unlist)
        return (GRanges())
      else
        return (GRangesList())
    }
      
    
    ff <- fread(paste("grep -v '^#'", x),sep='\t')
    if (ncol(ff)==10)
      setnames(ff, paste0("V",seq(1:10)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL"))
    else if (ncol(ff) == 11)
      setnames(ff, paste0("V",seq(1:11)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"))
    ff[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    ff$sample = gsub("(.*?)_.*","\\1",basename(x))
    ff[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    ff[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    ff[, NUMPARTS := as.integer(gsub(".*?NUMPARTS=([0-9]+).*", "\\1", INFO))]
    ff[, SCTG := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
    ff[, DISC_MAPQ := as.numeric(gsub(".*?DISC_MAPQ=([0-9]+).*", "\\1", INFO))]
    ff[, NM := as.integer(gsub(".*?;NM=([0-9]+).*", "\\1", INFO))]
    ff[, MATENM := as.integer(gsub(".*?;MATENM=([0-9]+).*", "\\1", INFO))]
    ff[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    print("...still formatting")
    ff[, REPSEQ := gsub(".*?;REPSEQ=([A-Z]+).*", "\\1", INFO)]
    ff[, REPSEQ := ifelse(grepl(";", REPSEQ), "", REPSEQ)] 
    ff[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
    ff[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    ff[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
    ff[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    if ("TUMOR" %in% colnames(ff)) {
      ff[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
      ff[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
      ff[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
    }
    if ("NORMAL" %in% colnames(ff)) {
      ff[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
      ff[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
      ff[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
    }
    ff[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
    ff[, inv := strand[1] == strand[2], by=uid]
    ff[, altstrand := rev(strand), by=uid]
    ff[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
    ff[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
    ff[, end := start]

    ff[, c("ID","REF","ALT","QUAL", "INFO", "GENO") := NULL]
    
    bad.ix <- ff[grepl("^G|^M", seqnames), uid]
    ff <- ff[!uid %in% bad.ix]
    
    ## if (length(aa) == 0) {
    ##   if (unlist)
    ##     return (GRanges())
    ##   else
    ##     return (GRangesList())
    ## }
    ## bad.ix <- unique(grl.unlist(aa)[grepl("^G", as.character(seqnames(grl.unlist(aa))))]$grl.ix)
    ## #tab <- table(grl.unlist(aa)$grl.ix[grepl("ref", as.character(seqnames(grl.unlist(aa))))])
    ## #bad.ix <- unique(c(bad.ix, as.numeric(names(tab[tab == 2]))))
    ## if (length(bad.ix)) 
    ##   aa <- aa[-bad.ix]
    ## if (length(aa))
    ##   mcols(aa)$filename = basename(x);
    ## if (unlist)
    ##   return (grl.unlist(aa))
    ## else
    ##   return (aa)

    df <- as.data.frame(dt2gr(ff))[,setdiff(colnames(ff), c("seqnames","start","end","strand"))]
    df <- df[!duplicated(df$uid),]
    g <- dt2gr(ff[,.(seqnames, start, end, strand, uid)])
    g <- split(g, g$uid)
    stopifnot(all(elementLengths(g) == 2))
    mcols(g) <- df[match(names(g), df$uid),]
    return(g)
  })
  return (snow4)
}

load_indel_dt <- function(x) {

  if (!file.exists(x))
    return(data.table())
  if (as.numeric(system(paste("grep -v '^#'", x, "| wc -l"), intern=TRUE)) == 0)
  return (data.table())

  suppressWarnings(ff <- fread(paste("grep -v '^#'", x),sep='\t'))

  if (!nrow(ff) || nrow(ff) > 20000)
    return(data.table())

  ff$start = ff$end = ff$V2
  setnames(ff, c("V1","V5","V10","V11","V4"), c("seqnames","ALT","normal","tumor","REF"))
  ff$sample = gsub("(.*?)_.*","\\1",basename(x))
  ff[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",V8))]
  ff[, INFO := gsub("READNAMES=(.*?);","",V8)]
  ff[, ID := paste(sample, V3,sep="_")]
  ff[, c("V8","V6","V9","V7","V2","V3") := NULL]
  ff[, del := nchar(REF) > nchar(ALT)]
  ff[, ins := nchar(ALT) > nchar(REF)]  
  ff[, SCTG := gsub(".*?SCTG=(.*?);.*","\\1",INFO)]
  ff[, REPSEQ := gsub(".*?REPSEQ=(.*?);.*","\\1",INFO)]
  ff[, DBSNP := grepl("DBSNP",INFO)]
  ff[, LOD := as.numeric(gsub(".*?LOD=(.*?);.*","\\1",INFO))]
  ff[, MAPQ := as.integer(gsub(".*?MAPQ=(.*?);.*","\\1",INFO))]
  ff[, NM := as.integer(gsub(".*?NM=(.*?);.*","\\1",INFO))]
  ff[, INFO := NULL]
  return(ff)
}

load_indel <- function(files, mc.cores=1) {

  snow4 <- parallel::mclapply(files, function(x) {
    print(paste(basename(x), "is", match(x, files), "of", length(files)))

    if (!file.exists(x)) {
      print(paste("File does not exist",x))
      return (GRanges())
    }
    if (exists('rowRanges'))
        fff <- rowRanges(rv <- VariantAnnotation::readVcf(x, "hg19"))
    else
        fff <- rowData(rv <- VariantAnnotation::readVcf(x, "hg19"))
    

    if (length(fff) == 0)
      return (GRanges())

    ### make sure we don't take germilne variants from strelka
    if ("NT" %in% colnames(info(rv)))
      non_norm_het = grepl("ref",info(rv)$NT)
    else
      non_norm_het <- rep(TRUE, length(fff))
    
    fff <- fff[ix <- (fff$FILTER == "PASS" | fff$FILTER == "alleleBias" | fff$FILTER == "QSI_ref" | fff$FILTER == "LOWLOD") & non_norm_het]
    if (length(fff) == 0)
      return (GRanges())
    
    mcols(fff)$filename = basename(x);

    mcols(fff)$REFWIDTH <- width(fff$REF)
    mcols(fff)$REF <- NULL
    tryCatch({mcols(fff)$ALTWIDTH <- unlist(nchar(mcols(fff)$ALT))}, error=function(e) { print ("ERROR in ALTWIDTH") }) #sapply(fff$ALT, nchar)
    mcols(fff)$ALT <- NULL
    
    if ("AD" %in% names(VariantAnnotation::geno(rv)) && length(dim(VariantAnnotation::geno(rv)$AD)) < 3 && ncol(VariantAnnotation::geno(rv)$AD) ==2) {
      if (is.list(VariantAnnotation::geno(rv)$AD[,2]))
        fff$AD <- sapply(VariantAnnotation::geno(rv)$AD[ix,2], function(x) x[2])
      else
        fff$AD <- VariantAnnotation::geno(rv)$AD[ix,2]
    }

    if ("AD" %in% names(VariantAnnotation::geno(rv)) && length(dim(VariantAnnotation::geno(rv)$AD)) < 3 && ncol(VariantAnnotation::geno(rv)$AD) == 1) {
      if (is.list(VariantAnnotation::geno(rv)$AD[,1]))
        fff$AD <- sapply(VariantAnnotation::geno(rv)$AD[ix,2], function(x) x[2])
      else
        fff$AD <- VariantAnnotation::geno(rv)$AD[ix,1]
    }

    
    if ("DP" %in% names(VariantAnnotation::geno(rv)) && length(dim(VariantAnnotation::geno(rv)$DP)) < 3 && !"QSI" %in% colnames(VariantAnnotation::info(rv)) && ncol(VariantAnnotation::geno(rv)$DP) == 2) {
      if (is.list(VariantAnnotation::geno(rv)$DP[,2]))
        fff$DP <- sapply(VariantAnnotation::geno(rv)$DP[ix,2], function(x) x[2])
      else
        fff$DP <- VariantAnnotation::geno(rv)$DP[ix,2]
    }
    
    if ("DP" %in% names(VariantAnnotation::geno(rv)) && length(dim(VariantAnnotation::geno(rv)$DP)) < 3 && !"QSI" %in% colnames(VariantAnnotation::info(rv)) && ncol(VariantAnnotation::geno(rv)$DP) == 1) {
      if (is.list(VariantAnnotation::geno(rv)$DP[,1]))
        fff$DP <- sapply(VariantAnnotation::geno(rv)$DP[ix,1], function(x) x[2])
      else
        fff$DP <- VariantAnnotation::geno(rv)$DP[ix,1]
    }

    ## platypus read counts
    if ("NV" %in% names(VariantAnnotation::geno(rv)) && length(dim(VariantAnnotation::geno(rv)$NV)) < 3) {
      nv_t <- sapply(VariantAnnotation::geno(rv)$NV[ix,2], max)
      nv_n <- sapply(VariantAnnotation::geno(rv)$NV[ix,1], max)
      names(nv_t) <- names(nv_n) <-  NULL
      #if (is.list(VariantAnnotation::geno(rv)$NV[,2]))
      fff$NV_n <- nv_n ##sapply(VariantAnnotation::geno(rv)$NV[ix,2], function(x) x[2])
      fff$NV_t <- nv_t ##sapply(VariantAnnotation::geno(rv)$NV[ix,2], function(x) x[2])
      fff$ratio <- ifelse(nv_n > 0, nv_t / nv_n, NA)  ##sapply(VariantAnnotation::geno(rv)$NV[ix,2], function(x) x[2])      
      #else
      #  fff$DP <- VariantAnnotation::geno(rv)$DP[ix,2]
    }

    ## strelka read counts
    print(names(VariantAnnotation::geno(rv)))
    if ("TIR" %in% names(VariantAnnotation::geno(rv))) {
      tt <- VariantAnnotation::geno(rv)$TIR[,"TUMOR",1]
      nn <- VariantAnnotation::geno(rv)$TIR[,"NORMAL",1]
      fff$TALT <- tt[ix]
      fff$NALT <- nn[ix]
    }

    
    if (sum(c("DP","AD") %in% colnames(mcols(fff)))==2) 
      fff$AF <- ifelse(fff$DP > 0, fff$AD/fff$DP, NA)

    #mcols(fff)$SPAN <- pmax(width(fff), width(unlist(fff$ALT)))
    mcols(fff)$CALC_SPAN <- pmax(width(fff), fff$ALTWIDTH - 1)
    mcols(fff)$ETYPE <- ifelse(fff$REFWIDTH > 1, 'del', 'ins')
    
    mcols(fff) <- cbind(VariantAnnotation::info(rv)[ix, ], mcols(fff))
    
    return (fff)
  }, mc.cores=mc.cores)
  return (snow4)
}

load_delly <- function(files) {

  delly4 <- lapply(files, function(x) {
    print(basename(x))
    if (!file.exists(x)) {
      print(paste("File does not exist",x))
      return (GRangesList())
    }
    dvcf <- readVcf(x, "hg19")

    if (exists("rowRanges"))
      d <- rowRanges(dvcf)
    else
      d <- rowData(dvcf)
      
    if (length(d) == 0)
      return (GRangesList())
    f = d$FILTER

    mcols(d) = info(dvcf)

    if (ncol(geno(dvcf)$DV)==2) {
      delly <- GRanges(c(as.character(seqnames(d)), as.character(d$CHR2)), IRanges(c(start(d), d$END), width=1), id=paste0("A",rep(seq_along(d), 2)), filter=rep(f,2),
                       SVTYPE=rep(mcols(d)$SVTYPE, 2), CT=rep(mcols(d)$CT, 2), TDISC=rep(geno(dvcf)$DV[,1],2),
                       NDISC=rep(geno(dvcf)$DV[,2],2),
                       TSPLIT=rep(geno(dvcf)$RV[,1],2), NSPLIT=rep(geno(dvcf)$RV[,2],2))
    } else {
      delly <- GRanges(c(as.character(seqnames(d)), as.character(d$CHR2)), IRanges(c(start(d), d$END), width=1), id=paste0("A",rep(seq_along(d), 2)), filter=rep(f,2),
                       SVTYPE=rep(mcols(d)$SVTYPE, 2), CT=rep(mcols(d)$CT, 2), DISC=rep(geno(dvcf)$DV[,1],2),
                       SPLIT=rep(geno(dvcf)$RV[,1],2))
    }
    
    #delly <- GRanges(c(as.character(seqnames(d)), as.character(d$MATECHROM)), IRanges(c(start(d), d$MATEPOS), width=1), id=paste0("A",rep(seq_along(d), 2)), filter=rep(f,2),
    #                 SVTYPE=rep(mcols(d)$SVTYPE, 2), CT=rep(mcols(d)$CT, 2))
    
    delly <- delly[delly$filter == "PASS"]

    if (ncol(geno(dvcf)$DV)==2) {
      dt <- data.table(id=delly$id, start=as.numeric(start(delly)), end=as.numeric(end(delly)), chr=as.character(seqnames(delly)),
                       SVTYPE=as.character(delly$SVTYPE), CT=as.character(delly$CT), TDISC=as.numeric(delly$TDISC), NDISC=as.numeric(delly$NDISC),
                       TSPLIT=as.numeric(delly$TSPLIT), NSPLIT=as.numeric(delly$NSPLIT))
    } else {
      dt <- data.table(id=delly$id, start=as.numeric(start(delly)), end=as.numeric(end(delly)), chr=as.character(seqnames(delly)),
                       SVTYPE=as.character(delly$SVTYPE), CT=as.character(delly$CT),
                       DISC=as.numeric(delly$DISC),
                       SPLIT=as.numeric(delly$SPLIT))
    }
    
    dt[, span := ifelse(chr[1] != chr[2], 1e9, abs(start[1] - end[2])), by=id]
    setkey(dt, id)
    dt <- unique(dt)
    
    del <- split(delly, delly$id)
    mcols(del)$id <- unique(grl.unlist(del)$id)
    mcols(del)$SPAN <- dt[mcols(del)$id]$span
    mcols(del)$SVTYPE <- dt[mcols(del)$id]$SVTYPE
    mcols(del)$CT <- dt[mcols(del)$id]$CT
    if (ncol(geno(dvcf)$DV)==2) {
      mcols(del)$TDISC <- dt[mcols(del)$id]$TDISC
      mcols(del)$NDISC <- dt[mcols(del)$id]$NDISC
      mcols(del)$TSPLIT <- dt[mcols(del)$id]$TSPLIT
      mcols(del)$NSPLIT <- dt[mcols(del)$id]$NSPLIT
    } else {
      mcols(del)$DISC <- dt[mcols(del)$id]$DISC
      mcols(del)$SPLIT <- dt[mcols(del)$id]$SPLIT
    }
    gg <- grl.unlist(del)
    ix <- mcols(gg)$CT == "3to5"
    strand(gg)[ix] <- rep(c("-", "+"), sum(ix)/2)
    ix <- mcols(gg)$CT == "5to3"
    strand(gg)[ix] <- rep(c("+", "-"), sum(ix)/2)
    ix <- mcols(gg)$CT == "3to3"
    strand(gg)[ix] <- rep(c("-", "-"), sum(ix)/2)
    ix <- mcols(gg)$CT == "5to5"
    strand(gg)[ix] <- rep(c("+", "+"), sum(ix)/2)
    gix <- gg$grl.ix
    mcols(gg) <- NULL
    gg <- split(gg, gix)
    
    ## assign the orientaitons
    ## ct <- mcols(del)$CT
    ## dd <- lapply(seq_along(del), function(i) {
    ##   b <- del[[i]]
    ##   if (ct[i] == "3to5") { ## del
    ##     strand(b) <- c('-', '+')
    ##   } else if (ct[i] == "5to3") { ## dup
    ##     strand(b) <- c('+', '-')
    ##   } else if (ct[i] == "3to3") {
    ##     strand(b) <- c("-", "-")
    ##   } else if (ct[i] == "5to5") {
    ##     strand(b) <- c("+", "+")        
    ##   }
    ##   return (b)
    ## })
    #dd <- GRangesList(dd)
    mcols(gg) <- mcols(del)
    return (gg)
})
  return (delly4)
}

rar_overlaps <- function(x, pad=500) {
  ro <- ra.overlaps(x, cons, pad=pad, ignore.strand=TRUE)
  fn <- cons[setdiff(seq_along(cons), ro[, 2])]
  fp <- x[setdiff(seq_along(x), ro[, 1])]
  paste(c("Num FN:", length(fn), "Num FP:", length(fp), "Num TP:", length(unique(ro[,2]))), collapse = " ")
}

indel_overlaps <- function(x, prefix="") {
  fo <- gr.findoverlaps(x + 10, gr.indels + 10,ignore.strand=TRUE)
  i_fn <- gr.indels[setdiff(seq_along(gr.indels), fo$subject.id)]
  i_fp <- x[setdiff(seq_along(x), fo$query.id)]
  paste(c(prefix,"Num FN:", length(i_fn), prefix,"Num FP:", length(i_fp), prefix,"Num TP:", length(unique(fo$query.id))))
}

flag.plot <- function(xindel = NULL, xsv = NULL, e, fname="plot.pdf", type="all", pad=10) {

  df <- data.frame()
  FPi <- FPs <- c()
  TPi <- TPs <- c()
  FNi <- FNs <- c()  

  e$span[e$span == -1] <- 1e8
  
  if (type == "medium") {
    num_events = length(unique(e$id[ix <- e$span >= 50 & e$span <= 200 & !grepl("RAR", e$id)]))
  } else if (type == "sv") {
    num_events = length(unique(e$id[ix <- e$span >= 500]))
  } else if (type == "indel") {
    num_events = length(unique(e$id[ix <- e$span < 50 & e$span > 0]))
  } else {
    num_events = length(unique(e$id[ix <- rep(TRUE, length(e))]))
  }
  
  if (!is.null(xindel)) {

    num_indel <- length(unique(e$id[e$span < 50 & e$span > 0]))
    er <- gr2dt(e)
    er[, subject.id := seq(nrow(er))]
    setkey(er, subject.id)
    
    suppressWarnings(fo <- gr2dt(gr.findoverlaps(xindel, e + 15)))
    setkey(fo, subject.id)

    fo <- er[,.(id, span, type, subject.id)][fo]
    fo[, qspan := xindel$SPAN[query.id]]
    fo[, diff := abs(span-qspan), by=id]
    fo <- fo[diff <= 10]
    
    TPi <- unique(fo$id)
    
    FNi <- setdiff(unique(e$id), unique(fo$id))
    FPi <- xindel[unique(setdiff(seq_along(xindel), fo$query.id))]
    ## only count false negatives among ones we should have seen
    FNi <- e[e$id %in% FNi & grepl("(del)|(ins)", e$id)]
    
    if (file.exists("/dev/shm/")) {
      saveRDS(FPi,"/dev/shm/fpi.rds")
      saveRDS(FNi,"/dev/shm/fni.rds")
    }
    
    #print(paste(c("TP Indel (<50bp)", length(TPi), "FP Indel", length(FPi))))
    #print(paste(c("Indel Precision: ", pr <- length(TPi)/(length(TPi) + length(FPi)), "Indel Recall:", rc<-length(TPi)/num_indel), collapse=" "))
    
  }

  if (!is.null(xsv)) {

    ## make GRangesList of true events
    grl.e <- split(e[ix], e$id[ix])
    this_id <- unique(grl.unlist(grl.e)$id)
    num_sv = length(unique(e$id[e$span >= 500]))
    
    ## get the overlaps
    id = seq_along(xsv)
    suppressWarnings(ro <- ra.overlaps(xsv, grl.e, pad=pad, ignore.strand=TRUE))

    FPs <- xsv[unique(setdiff(seq_along(xsv), id[ro[,1]]))]
    TPs <- unique(this_id[ro[,2]])
    FNs <- this_id[setdiff(seq_along(grl.e), ro[,2])]

    if (file.exists("/dev/shm/")) {
      saveRDS(FPs,"/dev/shm/fps.rds")
      saveRDS(FNs,"/dev/shm/fns.rds")
    }

  }

  print(paste("TP",sum(e$id %in% c(TPs,TPi) & e$span < 50)/2, sum(e$id %in% c(TPs,TPi) & e$span < 300 & e$span > 50)/2, sum(e$id %in% c(TPs,TPi) & e$span >= 500)/2))

  FPsprint1 <- tryCatch({sum(mcols(FPs)$SPAN < 50)}, error=function(e){0}) + ifelse(length(FPi), sum(FPi$SPAN < 50), 0)
  FPsprint2 <- tryCatch({sum(mcols(FPs)$SPAN > 50 & mcols(FPs)$SPAN < 300)}, error=function(e) { 0 }) + ifelse(length(FPi), sum(FPi$SPAN > 50 & FPi$SPAN < 300), 0)
  FPsprint3 <- tryCatch({sum(mcols(FPs)$SPAN >= 500)},error=function(e){0}) + ifelse(length(FPi), sum(FPi$SPAN >= 500), 0)
  
  print(paste("FP", FPsprint1, FPsprint2, FPsprint3))
  print(paste("  ",sum(e$span < 50), sum(e$span >= 50 & e$span <= 300), sum(e$span >= 500)))    
  #print(paste(c("--TP SV (>= 500)",  sum(e$id %in% c(TPs,TPi) & e$span >= 500)/2)))
  #print(paste(c("--TP Indel (< 50)", sum(e$id %in% c(TPs,TPi) & e$span < 50)/2)))
  #print(paste(c("--TP Med (50-300)", sum(e$id %in% c(TPs,TPi) & e$span < 300 & e$span > 50)/2)))
  #print(paste(c("--FP SV (>= 500)",  sum(mcols(FPs)$SPAN >= 500) + sum(FPi$SPAN >= 500))))
  #print(paste(c("--FP Indel (< 50)", sum(mcols(FPs)$SPAN < 50) + sum(FPi$SPAN < 50))))
  #print(paste(c("--FP Med (50-300)", sum(mcols(FPs)$SPAN > 50 & mcols(FPs)$SPAN < 300) + sum(FPi$SPAN > 50 & FPi$SPAN < 300))))

  FPs <- tryCatch({mcols(FPs)$SPAN},error=function(e){numeric(0)})
  dt <- data.table(TP=c(sum(e$id %in% c(TPs,TPi) & e$span >= 500)/2,
                        sum(e$id %in% c(TPs,TPi) & e$span < 50)/2,
                        sum(e$id %in% c(TPs,TPi) & e$span < 300 & e$span > 50)/2
                   ),
                   FP=c(sum(FPs >= 500) + sum(FPi$SPAN >= 500),
                        sum(FPs < 50) + sum(FPi$SPAN < 50),
                        sum(FPs > 50 & FPs < 300) + sum(FPi$SPAN > 50 & FPi$SPAN < 300)
                   ),
                   Group=c("SV","Indel","Medium")
  )
  
  tFP <- length(FPs) + length(FPi)
  tTP = length(unique(c(TPi,TPs)))
  #print(paste("Total FP:", tFP))
  #print(paste("Total TP:", tTP))
  #print(paste(c("Total Precision: ", pr <- tTP/(tTP + tFP), "Total Recall:", rc<-tTP/num_events), collapse=" "))
  
  if (fname=="noplot")
    return(dt)
  
  TPe <- e$id %in% c(TPs, TPi)
  FNe <- !TPe ##e$id %in% c(FNs, FNi)
  FP <- grbind(FPi, FPs)

  ## MAKE THE INDEL PLOT
  INDEL_CUTOFF = 9;
  ic.ix <- e$span == 1e8
  isindel <- e$span <= INDEL_CUTOFF 
  df.indel <- data.frame(span=c(e$span[TPe & isindel & !ic.ix], e$span[FNe & isindel & !ic.ix], FP$SPAN[FP$SPAN <= INDEL_CUTOFF & FP$SPAN > 0]),
                         type=c(rep("TP", sum(TPe & isindel & !ic.ix)), rep("FN", sum(FNe & isindel & !ic.ix)), rep("FP", length(FP$SPAN[FP$SPAN <= INDEL_CUTOFF & FP$SPAN > 0]))))
  df.indel$type = as.character(df.indel$type);    
  df.indel$type <- factor(df.indel$type, levels=c("TP", "FP", "FN"));   

  ## MAKE THE SV PLOT
  df.sv       <- data.frame(span=log10(c(e$span[TPe & !isindel], e$span[FNe & !isindel], FP$SPAN[FP$SPAN > INDEL_CUTOFF & FP$SPAN > 0])),
                            type=c(rep("TP", sum(TPe & !isindel)), rep("FN", sum(FNe & !isindel)), rep("FP", length(FP$SPAN[FP$SPAN > INDEL_CUTOFF & FP$SPAN > 0]))))
  df.sv$type = as.character(df.sv$type)
  df.sv$type <- factor(df.sv$type, levels=c("TP", "FP", "FN"))

   ## MAKE THE IC PLOT
   df.ic <- df.sv[df.sv$span==8,]
   df.ic$span <- factor(df.ic$span)

  labs <- c("TP"="True Positive", "FP"="False Positive", "FN"="False Negative")
  cols <- c("TP"="darkgreen", "FP"="darkred", "FN"="darkgrey")
  
  ## top, left, bottom, right
  g.sv <- ggplot() + geom_histogram(data=df.sv[df.sv$span < 7.7,],    aes(x=span, fill=type),binwidth=0.05) + theme_bw() + xlab("Distance (bp)") + ylab("Num Events") + scale_x_continuous(expand = c(0, 0), breaks=1:7, label=parse(text=paste("10", 1:7,sep="^"))) + scale_fill_manual(values=cols, name="", labels=labs) + coord_cartesian(xlim=c(1,7.7)) + theme(plot.margin=grid::unit(c(1,0,1,-0.5), "cm"), legend.position="none") + ylab("") + scale_y_continuous(expand = c(0, 0))
  g.in <- ggplot() + geom_histogram(data=df.indel, aes(x=span, fill=type),binwidth=1)   + theme_bw() + xlab("Distance (bp)") + ylab("Num Events") + scale_x_continuous(expand = c(0, 0), breaks=1:9) + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=cols, name="", labels=labs) + theme(legend.position="none",plot.margin=grid::unit(c(1,0.3,1,1), "cm")) + coord_cartesian(ylim=c(0,19600))
  g.ic <- ggplot() + geom_histogram(data=df.ic, aes(x=as.numeric(span), fill=type), binwidth=1) + theme_bw() + xlab("") + ylab("Num Events") + scale_x_discrete(label=c("IC"),expand = c(0, 0)) + theme(plot.margin=grid::unit(c(1,0.3,1.5,-0.3), "cm"),legend.position="none") + ylab("") + scale_fill_manual(values=cols, name="", labels=labs) + scale_y_continuous(expand = c(0, 0))
  require(gridExtra)

  g <- grid.arrange(g.in, g.sv, g.ic, ncol=3, widths=c(2.5,5,0.8)) #, width=8, height=3, filename=fname)

  return(list(df=df, g=g, TP=TPe, FPi=FPi, FPs=FPs))
}

grl2links <- function(x) {

  y = grl.unlist(x);

  ix = y$grl.iix == 1;
  sn = as.character(seqnames(y))
  p = as.numeric(start(y));
  links = data.frame(Chromosome=as.character(sn[ix]), chromStart=p[ix], chromEnd=p[ix], Chromsome.1=as.character(sn[!ix]), chromStart.1=p[!ix], chromEnd.1=p[!ix], stringsAsFactors=FALSE)
  return (links)
  
}

grl2circos <- function(x, genes = GRanges(), chr.exclude=NULL) {

  library(RCircos)
  data(UCSC.HG19.Human.CytoBandIdeogram);
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
  tracks.inside <- 10;
  tracks.outside <- 0;
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);

  ## get the gene label dat
  gene.dat <- data.frame()
  if (length(genes)) {
    genes <- genes[width(genes) < 2e6]
    fo1 <- gr.findoverlaps(gr1+10e3, genes)
    fo2 <- gr.findoverlaps(gr2+10e3, genes)
    
    fo <- c(fo1,fo2)
    gene.dat <- data.frame()
    if (length(fo)) {
      fo <- fo[!duplicated(fo$subject.id)]
      gene.dat = data.frame(Chromosome = seqnames(genes[fo$subject.id]), chromStart=start(genes[fo$subject.id]),
        chromEnd=end(genes[fo$subject.id]), Gene=genes$gene[fo$subject.id])
      print(gene.dat)
    }
  }
  gename.col <- 4;
  side <- "in";
  track.num <- 1;

  ## set the links data
  links <- grl2links(x)

  ## draw
  RCircos.Set.Plot.Area();
  RCircos.Chromosome.Ideogram.Plot();
  if (nrow(gene.dat) > 0) {
    RCircos.Gene.Connector.Plot(gene.dat, track.num, side);
    track.num <- 2;
    name.col <- 4;
    RCircos.Gene.Name.Plot(gene.dat, name.col,track.num, side);
  }
  
  if (nrow(links) > 0)
    RCircos.Link.Plot(links, track.num, by.chromosome=TRUE) ## by.chromosome is for color

}

evaluate_complex <- function() {


}

load_pindel <- function(x, tum = 6, normalt =0, normref = 10) {

  ff <- fread(x, sep='\t')
  ff[, END := as.numeric(gsub("END=([0-9]+);.*", "\\1", V4))]
  #ff[, TUM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V6))]
  #ff[, TUM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V6))]
  #ff[, NORM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V7))]
  #ff[, NORM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V7))]
  ff[, HOMSEQ := ifelse(grepl("HOMSEQ", V4), nchar(gsub("HOMSEQ=([A-Z]+);.*", "\\1", V4)), 0)]
  ff[, TUM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V7))]  ## this is true for HCC1143_mem
  ff[, TUM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V7))]
  ff[, NORM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V6))]
  ff[, NORM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V6))]
  
  ff[, SPAN := abs(as.numeric(gsub(".*?SVLEN=(-?[0-9]+).*","\\1",ff$V4)))]
  pindel10x_s <- with(ff[NORM_ALT <= normalt & TUM_ALT >= tum & NORM_REF >= normref & HOMSEQ < 30 & TUM_ALT / TUM_REF >= 0.2], GRanges(V1, IRanges(as.numeric(V2),END), SPAN=SPAN, END=END, TUMALT=TUM_ALT, NORMAL=NORM_ALT, TUMREF=TUM_REF, NORMREF=NORM_REF, FILTER=V3, HOMSEQ=HOMSEQ))
  pindel10x_g <- with(ff[NORM_ALT > 2 | TUM_ALT / NORM_ALT < 10], GRanges(V1, IRanges(as.numeric(V2),END), SPAN=SPAN))
  
  fo <- gr.findoverlaps(pindel10x_g + 2, pindel10x_s)
  #pindel10x_s2 <- pindel10x_s[setdiff(seq_along(pindel10x_s), fo$subject.id)]

  return (pindel10x_s)
  
}

load_pindel_vcf <- function(x, tum = 6, normalt =0, normref = 2, af=0.05) {

  ff <- fread(paste("grep -v '^#'", x),sep='\t')  
  ff[, END := as.numeric(gsub("END=([0-9]+);.*", "\\1", V8))]
  #ff[, TUM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V6))]
  #ff[, TUM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V6))]
  #ff[, NORM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V7))]
  #ff[, NORM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V7))]
  ff[, HOMSEQ := ifelse(grepl("HOMSEQ", V4), nchar(gsub("HOMSEQ=([A-Z]+);.*", "\\1", V8)), 0)]
  ff[, TUM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V11))]  ## this is true for HCC1143_mem
  ff[, TUM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V11))]
  ff[, NORM_ALT := as.numeric(gsub(".*?:[0-9]+,([0-9]+)", "\\1", V10))]
  ff[, NORM_REF := as.numeric(gsub(".*?:([0-9]+),[0-9]+", "\\1", V10))]
  
  ff[, SPAN := abs(as.numeric(gsub(".*?SVLEN=(-?[0-9]+).*","\\1",V8)))]
  pindel_s <- with(ff[NORM_ALT <= normalt & TUM_ALT >= tum & NORM_REF >= normref & HOMSEQ < 30 & TUM_ALT / TUM_REF >= af], GRanges(V1, IRanges(as.numeric(V2),END), SPAN=SPAN, END=END, TUMALT=TUM_ALT, NORMAL=NORM_ALT, TUMREF=TUM_REF, NORMREF=NORM_REF, FILTER=V3, HOMSEQ=HOMSEQ))
  #pindel10x_g <- with(ff[NORM_ALT > 2 | TUM_ALT / NORM_ALT < 10], GRanges(V1, IRanges(as.numeric(V2),END), SPAN=SPAN))
  
  #fo <- gr.findoverlaps(pindel10x_g + 2, pindel10x_s)
  #pindel10x_s2 <- pindel10x_s[setdiff(seq_along(pindel10x_s), fo$subject.id)]

  return (pindel_s)
  
}


load_freebayes <- function(x, tum = 4, normalt = 0, normref = 10) {

  ff <- fread(paste("grep -v '^#'", x),sep='\t')
  ff <- ff[!grepl("snp", ff$V8)]
  ff[, SSC := as.numeric(gsub(".*?SSC=(.*?);.*","\\1",V8))]
  #ff[, TUM_ALT := as.numeric(gsub("[0-9]+/[0-9]+:.*?:.*?:.*?:.*?:([0-9]+):.*", "\\1", V11))]
  #ff[, TUM_REF := as.numeric(gsub(".*?:.*?:.*?:([0-9]+).*", "\\1", V11))]
  #ff[, NORM_ALT := as.numeric(gsub(".*?:.*?:.*?:.*?:.*?:([0-9]+).*", "\\1", V10))]
  #ff[, NORM_REF := as.numeric(gsub(".*?:.*?:.*?:([0-9]+).*", "\\1", V10))]
  ff[, SPAN := abs(as.numeric(gsub(".*?LEN=([0-9]+).*","\\1",ff$V8)))]
  #out <- with(f2 <- ff[ff$NORM_ALT <= normalt & ff$TUM_ALT >= tum & ff$NORM_REF >= normref], GRanges(V1, IRanges(as.numeric(V2),as.numeric(V2)), SPAN=SPAN))
  #out <- with(f2 <- ff[ff$SSC > 25], GRanges(V1, IRanges(as.numeric(V2),as.numeric(V2)), SPAN=SPAN))  

  ## GL that its ref
  ff[, REF_TUMOR_GL := {
    ab <- lapply(strsplit(V11,':'), function(x) strsplit(x[8], ",")[[1]]);
    return(ifelse(sapply(ab, length)==1, -1e6,-sapply(ab, function(y) {
      return(as.numeric(y[1]))
    })))
  }]

  ## parse the genotype likelihoods as provided by FreeBayes
  ff[, MAX_NONREF_TUMOR_GL := {
    ab <- lapply(strsplit(V11,':'), function(x) strsplit(x[8], ",")[[1]]);
    return(ifelse(sapply(ab, length)==1, -1e6, sapply(ab, function(y) {
      tt <- as.numeric(y[1]);
      if (length(y) > 1)
        a2 <- max(sapply(y[2:length(y)], function(z) as.numeric(z) - tt))
      else
        return(-1)
      return(a2)
    })))
  }]

  ## parse the genotype likelihoods as provided by FreeBayes
  ff[, MAX_NORMAL_GL := {
    ab <- lapply(strsplit(V10,':'), function(x) strsplit(x[8], ",")[[1]]);
    return(unlist(ifelse(sapply(ab, length) == 1, -1e6, sapply(ab, function(y) {
      tt <- as.numeric(y[1]);
      if (length(y) > 1)
        return(min(sapply(y[2:length(y)], function(z) tt - as.numeric(z))))
      return(-1)
    }))))
  }]

  ff[, TUMALT := as.numeric(gsub("[0-9]/[0-9]:.*?:.*?:.*?:.*?:(.*?):.*", "\\1", V11))]
  ff[, NORMALT := as.numeric(gsub("[0-9]/[0-9]:.*?:.*?:.*?:.*?:(.*?):.*", "\\1", V10))]  

  out <- with(ff[MAX_NONREF_TUMOR_GL > 8 & MAX_NORMAL_GL > 8 & NORMALT <= 2 & TUMALT >= 4 & ifelse(NORMALT > 0, TUMALT / NORMALT, 100) > 20], GRanges(V1, IRanges(as.numeric(V2),as.numeric(V2)), SPAN=SPAN))  
  
}


span.hist <- function(s) {

  s[s < 0] <- 1e9
  s <- s[s >= 10]
  s[ix <- s == 1e9] <- 1e9 + runif(sum(ix)) * 10e9
  #df <- with(sanger, data.frame(s=log10(s), chr=factor(V1, levels=c(seq(22),"X")), type=paste(paste0(strand,altstrand), ifelse(sanger$intra_tad,"TAD","INTER-TAD")))
  df <- with(sanger[!duplicated(urar_id)], data.frame(s=log10(ifelse(SPAN!=-1, SPAN, 1e9)), chr=factor(V1, levels=c(seq(22),"X")), type=paste0(strand,altstrand), per_sample_count=per_sample_count))
  df <- data.frame(s=log10(s), chr=factor(gr.snow$seqnames, levels=c(seq(22),"X")), type=paste0(gr.snow$strand,gr.snow$altstrand), valid=gr.snow$valid & !gr.snow$germbad & !(gr.snow$in_line & gr.snow$invtype))
  
  ## tmp
  ## df <- snow_B_mini
  ## df$span <- df$SPAN
  ## df$span[df$span < 0] <- 1e9
  ## df <- df[df$span >= 10]
  ## df[, s := log10(span)]

  df <- df[gr.sanger$tad != 0 & gr.sanger$SPAN < 1e9,]
  df <- df[df$valid,]
  df$event_count <- ""
  df$event_count[df$per_sample_count < 20] = "<10"
  df$event_count[df$per_sample_count >= 11 & df$per_sample_count <= 25] = "[11,25]"
  df$event_count[df$per_sample_count >= 26 & df$per_sample_count <= 50] = "[26,50]"
  df$event_count[df$per_sample_count >= 51 & df$per_sample_count <= 100] = "[51,100]"
  df$event_count[df$per_sample_count >= 101 & df$per_sample_count <= 300] = "[101,300]"
  df$event_count[df$per_sample_count >= 301 & df$per_sample_count <= 10000] = "[301,10000]"
  df$event_count[df$per_sample_count > 10000] = ">10000"          
  #df$type <- "NON_INV"
  #df$type[df$inv] <- "INV"

  g <- ggplot() + geom_histogram(data=df, aes(x=s, fill=type), alpha=0.5, binwidth=0.1, position='dodge') +theme_bw() + xlab("Distance (bp)") + ylab("Num Events") + scale_x_continuous(breaks=1:8, label=parse(text=paste("10", 1:8,sep="^"))) +
    facet_wrap(~ event_count, scale='free') + scale_fill_manual(values=c("darkgreen", "blue", "red", "green"))
                                         #geom_histogram(data=data.frame(width=width(tad)), aes(x=log10(width)), binwidth=0.1)
                                        #+ facet_wrap(~ chr)##+ theme(legend.position="none")

  skitools::ppdf(print(g), width=12, height=12, "sanger_hist_bycount.pdf")
}

bar.plot <- function(s) {

  tab <- table(s)
  tab <- tab[order(tab, decreasing=TRUE)]
  df <- data.frame(names=factor(names(tab), levels=names(tab)), count=tab)
  g <- ggplot() + geom_bar(data=df, aes(x=names, y=count), stat='identity')
  return(g)
}

dump.to.bed <- function(gr, file='foo.bed') {

  options(scipen=999)
  #gr <- sort(gr)
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   #names=c(rep(".", length(gr))),
                   scores=gr$score,
                   #scores=c(rep(".", length(gr))),
                   strands=strand(gr))
  df <- cbind(df, mcols(gr))
  
  write.table(df, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
}

distance.to.nearest.plot <- function(ss) {

  ## distance between events
  p1 <- ss[!duplicated(urar_id)]$pos1
  c1 <- ss[!duplicated(urar_id)]$seqnames
  dist <- unlist(sapply(unique(c1), function(x) diff(sort(p1[c1==x]))))
  
  ## make random positions
  .getrands <- function(num) {
    rand_starts <- sample(sum(as.numeric(seqlengths(si)[1:23])), num)
    csum <- cumsum(as.numeric(seqlengths(si)[1:23]))
    chr <- sapply(rand_starts, function(x) { which(csum - x > 0)[1] })
    csum <- c(0,csum)
    real_starts <- sapply(seq_along(rand_starts), function(x) rand_starts[x] - csum[chr[x]])
    return(data.table(chr=chr, real_starts=real_starts))
  }

  grvalid <- gr.100map[gr.100map$score < 0.75]
  r1 <- .getrands(length(dist)*1.5) ## get more than needed, bc we trim later
  gr <- GRanges(r1$chr, IRanges(r1$real_starts, width=1))
  fo <- gr.findoverlaps(gr, grvalid)
  r1 <- r1[!fo$query.id]
  r1 <- r1[sample(nrow(r1), length(dist))] ## get it back to right size
  rand_distances <- unlist(sapply(unique(r1$chr), function(x) diff(sort(r1$real_starts[r1$chr==x]))))
  
  df <- data.frame(x=c(rand_distances, dist), group=c(rep("Random", length(rand_distances)), rep("Real", length(dist))))
  #g <- ggplot(df) + geom_freqpoly(aes(x=log10(x), color=group), binwidth=0.1, position='dodge') + xlab("Distance (bp)") + ylab("Num Events") +
  #  scale_x_continuous(breaks=1:8, label=parse(text=paste("10", 1:8,sep="^"))) + scale_color_manual(values=c("gray50", "darkblue"), name="") +
  #    theme(legend.position=c(0.8,0.5))
  #return(g)
  return(df)

}

scatter_plot <- function(xdat, ydat, xlabr, ylabr) {

  ix <- !is.na(xdat) & !is.na(ydat)
  xdat <- xdat[ix]
  ydat <- ydat[ix]  
  
  lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(coef(m)[1], digits = 2),
                          b = format(coef(m)[2], digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }

  df = data.frame(x=xdat, y=ydat);

  XSPOT = median(xdat)
  YSPOT = median(ydat)
  p <- ggplot(data=df,aes(x, y))+
    #stat_summary(fun.data=mean_cl_normal) +
    #geom_smooth(method='lm',formula=y~x) +
      geom_point() + coord_cartesian(xlim=c(0,max(xdat)), ylim=c(0,max(ydat)))
  p <- p + ggtitle(parse(text=lm_eqn(df))) #geom_text(x = XSPOT, y = YSPOT, label = lm_eqn(df), parse = TRUE) + xlab(xlabr) + ylab(ylabr)
  return(p)
  
}

plot.per.elem <- function(data, gr, val, by = NULL, NUM=1:20, scale.width=FALSE) {

  fo <- gr.findoverlaps(gr, data)
  fo$val <- mcols(gr)[,val][fo$query.id]
  if (!is.null(by)) {
    fo$by <- data[,.(eval(parse(text=by)))][[1]][fo$subject.id]
    fo <- fo[!duplicated(paste(fo$val, fo$by))]
  }

  tab <- table(fo$val)

  if (scale.width) {
    gs <- structure(names=mcols(gr)[,val], width(gr))
    sgs <- tab/gs[names(tab)] * 1e3 ## print by kb
  } else {
    sgs <- tab
  }
  
  ord <- rev(order(sgs))
  df <- data.frame(val=names(tab)[ord], count=as.numeric(tab)[ord], scale_count=sgs[ord]) 
  df$val <- factor(df$val, levels = df$val)
  
  g <- ggplot(data=df[NUM,]) + geom_bar(aes(x=val, y=scale_count), stat='identity', fill='gray80',color='black') + theme_bw() +
    xlab("") + theme(axis.text.x=element_text(angle=90))
      #ylab("Breakpoints per gene") + 
  #ppdf(print(g), width=7, height=4, "merged_event_count_per_gene.pdf")         
  
}

subset_rars <- function(gr, xxx, pad=0) {

  if (is.character(gr)) {
    if (gr %in% gr.genes$gene)
      gr <- gr.genes[gr.genes$gene == gr]
    else if (grepl(":", gr))
      gr <- gUtils::parse.grl(gr)[[1]]
  }

  gr <- gr + pad
  
  gg <- split(gg<-dt2gr(xxx[uid %in% unique((dt2gr(xxx) %*% gr)$uid)]), gg$uid)
  mc <- mcols(grl.unlist(gg))
  mc <- mc[mc$grl.iix==1,]
  mcols(gg) <- mc
  return(gg)
}

code_enrichment_score <- function(codes, sifs) {

  tab <- table(sifs$dcc_project_code)

  dt <- data.table(codes=codes, sifs.count=tab[codes])
  dt[, code.prob := sifs.count / nrow(sifs)]
  dt[, code.count := nrow(.SD), by=codes]
  setkey(dt, codes)
  dt <- unique(dt)
  total.codes = length(codes)

  dt[, enrichment.pval := pbinom(code.count, total.codes, p=code.prob, lower.tail=FALSE)]
  dt[, log10p := log10(enrichment.pval)]
  dt[, neglog10p := -log10p]
  setkey(dt, log10p)
  return(dt)
}

#pp <- fread("grep IN_BK /broad/broadsv/NA12878/PacBio/PacBioSplit.txt")
#setnames(pp, c("V2","V3", "V4"), c("seqnames", "start", "strand"))
#pp[, end := start]
#grp <- dt2gr(pp)
#mcols(grp) <- NULL
annotate_pacbio <- function(x, pad=10) {

  if (class(x) == "GRangesList") {
    gr <- grl.unlist(x)
    fo <- gr2dt(gr.findoverlaps(gr + pad, grp, max.chunk=1e10, mc.cores=5))
    fo[, qname := pp$V5[subject.id]]
    fo[, qname.count.per.bkp.id := sum(!duplicated(qname)), by=query.id]
    fo[, rar.id := gr$grl.ix[query.id], query.id]
    fo[, grl.iix := gr$grl.iix[query.id], query.id]
    fo[, qname.count.per.rar.id := sum(!duplicated(qname)), by=rar.id]
    fo[, grl.iix.count.per.rar.id := sum(!duplicated(grl.iix)), by=rar.id]
    r <- fo[qname.count.per.rar.id >= 2 & grl.iix.count.per.rar.id == 2]
    setkey(r, rar.id)
    mcols(x)$pacbio = FALSE
    mcols(x)$pacbio[unique(r)$rar.id] <- TRUE
    return(x)
  } else if (class(x) == "GRanges") {
    fo <- gr2dt(gr.findoverlaps(x + pad, grp, max.chunk=1e10, mc.cores=5))
    fo[, qname := pp$V5[subject.id]]
    fo[, qname.count.per.bkp.id := sum(!duplicated(qname)), by=query.id]
    fo[, rar.id := gr$grl.ix[query.id], query.id]
    fo[, grl.iix := gr$grl.iix[query.id], query.id]
    fo[, qname.count.per.rar.id := sum(!duplicated(qname)), by=rar.id]
    fo[, grl.iix.count.per.rar.id := sum(!duplicated(grl.iix)), by=rar.id]
    r <- fo[qname.count.per.rar.id >= 2 & grl.iix.count.per.rar.id == 2]
    setkey(r, rar.id)
    mcols(x)$pacbio = FALSE
    mcols(x)$pacbio[unique(r)$rar.id] <- TRUE
  }
  
}


enrichment_test <- function(data, track) {

  ##genome_size = sum(as.numeric(width(setdiff(setdiff(gr.stripstrand(si2gr(gUtils::si)), x@misc$blacklist), x@misc$hengli))))
  genome_size <- 2429413590

  fo <- gr.findoverlaps(data, track)
  if ("uid" %in% colnames(data)) {
    event_count <- length(unique(data$uid))
    ac <- length(unique(data$uid[fo$query.id]))
  } else {
    event_count = nrow(data)
    ac = length(unique(fo$query.id))
  }

  a <- ac / event_count
  bc <- sum(as.numeric(width(reduce(track))))
  b <- bc / genome_size
  enrich = a / b
  data.frame(min=  qbinom(0.025, event_count, a) / event_count / b, max=qbinom(0.975, event_count, a) / event_count / b, val=ac / event_count / b, enrich = enrich)
}

rep_extract <- function(x, elem) {

  return(x@misc$repmask[x@misc$repmask$repeat.element == elem])
  
}

enr_over_rand <- function(sites, rands, events, track = NULL) {

  ## sites .e.g. CTCF
  ## rands e.g. random CTCF
  ## events e.g. SV dels
  ## track e.g. tad boundaries
  N=length(rands)/length(sites)

  ## get sites covering track (e.g boundary CTCF)
  if (!is.null(track)) {
    print("...getting sites (e.g. CTCF) covering track (e.g. TAD boundaries)")
    fo <- gr.findoverlaps(sites, track)
    sites.ctad <- sites[unique(fo$query.id)]
  } else {
    sites.ctad = sites
  }

  ## find number of sites covering events
  print("...getting sites covering events (e.g. dels)")
  fo <- gr.findoverlaps(sites.ctad, events)
  real <- length(unique(fo$query.id))

  ## find rand sites covering track
  if (!is.null(track)) {
    print("...getting rand sites covering tracks (e.g. TAD boundaries)")
    fo <- gr.findoverlaps(track, rands, max.chunk=Inf, mc.cores=1)
    rands.ctad <- rands[unique(fo$subject.id)]
  } else {
    rands.ctad <- rands
  }
  
  ## find rand sites covering events
  print("...getting rand sites covering events")
  fo.r <- gr.findoverlaps(rands.ctad, events, max.chunk=Inf, mc.cores=1)
  fo.r$run = fo.r$query.id %% N

  ## count rand sites covering events (get distribution)
  print("...tabulating rand overlaps")
  tab <- table(fo.r$run[unique(fo.r$query.id)])

  ## return the plot
  mean=mean(as.numeric(tab))
  sd=sd(as.numeric(tab))
  return(list(
              g=ggplot() + geom_histogram(data=data.frame(x=as.numeric(tab) / (length(rands.ctad)/1000/length(sites.ctad))), aes(x=x), binwidth=10) + theme_bw() +
              ylab("Count") + geom_line(data=data.frame(x=c(real,real), y=c(0,50)), aes(x=x,y=y), color='red'),
              OR=real/mean(as.numeric(tab)),
              ORL=real/(mean + 1.96*sd), ORH=real / (mean - 1.96*sd)
              ))

}

span.histogram <- function(span) {

  delly <- mcols(load_delly("/broad/hptmp/jwala/delly/rearrangements.somatic.vcf")[[1]])$SPAN
  #delly <- mcols(load_delly(ind_c2$delly_vcf[1])[[1]])$SPAN  
  strelka <- mcols(load_indel(ind_c3$strelka_passed_somatic_indel_vcf_file_wgs[2])[[1]])$CALC_SPAN
  #strelka <- mcols(load_indel(ind_c2$strelka_passed_somatic_indel_vcf_file_wgs[1])[[1]])$CALC_SPAN  
  snowman <- mcols(load_snowman(ind_c3$snowman_somatic_vcf[2])[[1]])$SPAN
  #nowman <- mcols(load_snowman(ind_c2$snowman_somatic_vcf[1])[[1]])$SPAN  
  snowmani <- mcols(load_indel(ind_c3$snowman_somatic_indel_vcf[2])[[1]])$SPAN
  #snowmani <- mcols(load_indel(ind_c2$snowman_somatic_indel_vcf[1])[[1]])$SPAN  
  spans <- c(snowman, snowmani) #[sample(length(snowmani), length(snowmani)*0.7)])
  spand <- c(delly, strelka) #[sample(length(strelka), length(strelka)*0.7)])
  #spand <- c(spand, rep(1e9, length(spans) - length(spand)))
  df <- data.frame(span=c(spans, spand, gr.events$span[ix2 <- !duplicated(gr.events$id)]), group=c(rep("s", length(spans)),rep("d", length(spand)), rep("r", sum(ix2))))
  df$span[ix <- df$span < 10] <- df$span[ix] + runif(sum(ix))
  df$x <- log10(df$span)
  df$x[df$span == 1] <- 0.11
  df$x[df$span == 2] <- 0.21
  df$x[df$span == 3] <- 0.31
  df$x[df$span == 4] <- 0.41
  df$x[df$span == 5] <- 0.51    
  df$x[df$span == 6] <- 0.61
  df$x[df$span == 7] <- 0.71
  df$x[df$span == 8] <- 0.81
  df$x[df$span == 9] <- 0.91  
  
  g <- ggplot(df[!is.na(df$x) & df$span >= 1 & df$group %in% c('s') & runif(nrow(df)) < 0.5, ]) + geom_histogram(aes(x=x), fill='darkgreen', position='dodge', binwidth=0.1) + theme_bw() + #facet_wrap(~ group, nrow=2) +
  #g <- ggplot(df[!is.na(df$x) & df$span >= 3, ]) + geom_freqpoly(aes(x=x, color=group), alpha=1, adjust=1/3)+ theme_bw() + #facet_wrap(~ group, nrow=2) +     
    xlab("Distance (bp)") + ylab("Num Events") +
      scale_x_continuous(breaks=0:8, label=parse(text=paste("10", 0:8,sep="^"))) + coord_cartesian(xlim=c(0,8)) #,ylim=c(0,600))
  ppdf(print(g), width=5, height=3)  
  
}
