#!/usr/bin/env Rscript

## set the right library paths
.libPaths = c("/xchip/gistic/Jeremiah/R", "/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/lib64/R/library")

require(data.table)
ra_breaks <-
function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE,  breakpointer = FALSE, seqlevels = NULL,
    get.loose = FALSE ## if TRUE will return a list with fields $junctions and $loose.ends
  )
  {
      if (is.character(rafile))
          {
              if (grepl('(vcf$)|(vcf.gz$)', rafile))
                  {
                      library(VariantAnnotation)
                      vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))
                      if (!('SVTYPE' %in% names(info(vcf)))) {
                        warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                          return(GRangesList());
                        }
                      
                      vgr = rowData(vcf) ## parse BND format

                      ## no events
                      if (length(vgr) == 0)
                        return (GRangesList())

                      ## fix mateids if not included
                      if (!"MATEID"%in%colnames(mcols(vgr))) {
                        nm <- vgr$MATEID <- names(vgr)
                        ix <- grepl("1$",nm)
                        vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                        vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                        vgr$SVTYPE="BND"
                      }
                      
                      if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                        stop("MATEID or SVTYPE not included. Required")
                      vgr$mateid = info(vcf)$MATEID
                      vgr$svtype = info(vcf)$SVTYPE

                      if (!is.null(info(vcf)$SCTG))
                          vgr$SCTG = info(vcf)$SCTG
                      
                      if (sum(vgr$svtype == 'BND')==0)
                          stop('Vcf not in proper format.  Will only process rearrangements in BND format')

                      if (!all(vgr$svtype == 'BND'))
                          warning(sprintf('%s rows of vcf do not have svtype BND, ignoring these'), sum(vgr$svtype != 'BND'))

                      bix = which(vgr$svtype == "BND")
                      vgr = vgr[bix]
                      vgr$first = !grepl('^(\\]|\\[)', vgr$ALT) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                      vgr$right = grepl('\\[', vgr$ALT) ## ? are the (sharp ends) of the brackets facing right or left
                      vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                      vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', vgr$ALT))
                      vgr$mcoord = gsub('chr', '', vgr$mcoord)

                      if (all(is.na(vgr$mateid)))
                          if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                              {
                                  warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                                  vgr$mateid = paste(gsub('::\\d$', '', names(vgr)), (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                              }
                          else if (!is.null(vgr$SCTG))
                              {
                                  warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                                  require(igraph)
                                  ucoord = unique(c(vgr$coord, vgr$mcoord))
                                  vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                                  if (any(duplicated(vgr$mateid)))
                                      {
                                          warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                                          bix = bix[!duplicated(vgr$mateid)]
                                          vgr = vgr[!duplicated(vgr$mateid)]
                                      }
                              }
                          else
                              stop('MATEID tag missing')

                      vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                      pix = which(!is.na(vgr$mix))

                      vgr.pair = vgr[pix]

                      if (length(vgr.pair)==0)
                          stop('No mates found despite nonzero number of BND rows in VCF')
                      vgr.pair$mix = match(vgr.pair$mix, pix)
                      vix = which(1:length(vgr.pair)<vgr.pair$mix )
                      vgr.pair1 = vgr.pair[vix]
                      vgr.pair2 = vgr.pair[vgr.pair1$mix]

                      ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint
                      
                      ## if "first" and "right" then we set this entry "-" and the second entry "+"
                      tmpix = vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      ## if "first" and "left" then "-", "-"
                      tmpix = vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "left" then "+", "-"
                      tmpix = !vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "right" then "+", "+"
                      tmpix = !vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
                      if (any(pos1))
                          {
                              start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                              end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                          }

                      pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
                      if (any(pos2))
                          {
                              start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                              end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
                          }
                      ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

                      this.inf = info(vcf)[bix[pix[vix]], ]

                      if (is.null(this.inf$POS))
                          this.inf = cbind(data.frame(POS = ''), this.inf)
                      if (is.null(this.inf$CHROM))
                          this.inf = cbind(data.frame(CHROM = ''), this.inf)

                      if (is.null(this.inf$MATL))
                          this.inf = cbind(data.frame(MALT = ''), this.inf)
                      
                      this.inf$CHROM = seqnames(vgr.pair1)
                      this.inf$POS = start(vgr.pair1)
                      this.inf$MATECHROM = seqnames(vgr.pair2)
                      this.inf$MATEPOS = start(vgr.pair2)
                      this.inf$MALT = vgr.pair2$ALT
                      
                      values(ra) = cbind(fixed(vcf)[bix[pix[vix]],], this.inf)
                      
                      if (is.null(values(ra)$TIER))
                          values(ra)$tier = ifelse(values(ra)$FILTER == "PASS", 2, 3) ## baseline tiering of PASS vs non PASS variants
                      else
                          values(ra)$tier = values(ra)$TIER

                      if (!get.loose)
                          return(ra)
                      else
                          {
                              npix = is.na(vgr$mix)
                              vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation
                              values(vgr.loose) = cbind(fixed(vcf)[bix[npix], ], info(vcf)[bix[npix], ])

                              return(list(junctions = ra, loose.ends = vgr.loose))
                          }
                  }
              else
                  rafile = read.delim(rafile)
          }
            
     if (is.data.table(rafile))
         rafile = as.data.frame(rafile)

    if (nrow(rafile)==0)
        {
            out = GRangesList()
            values(out) = rafile
            return(out)
        }
  
    if (snowman) ## flip breaks so that they are pointing away from junction
      {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
      }
      
    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
      {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]        
      }

     
    if (is.null(rafile$str1))
      rafile$str1 = rafile$strand1

    if (is.null(rafile$str2))
      rafile$str2 = rafile$strand2
     if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
         {
             if (breakpointer)
                 {
                     rafile$pos1 = rafile$T_BPpos1
                     rafile$pos2 = rafile$T_BPpos2
                 }
             
             if (!is.numeric(rafile$pos1))
                 rafile$pos1 = as.numeric(rafile$pos1)

             if (!is.numeric(rafile$pos2))
                 rafile$pos2 = as.numeric(rafile$pos2)

             ## clean the parenthesis from the string

             rafile$str1 <- gsub('[()]', '', rafile$str1)
             rafile$str2 <- gsub('[()]', '', rafile$str2)

             if (is.character(rafile$str1) | is.factor(rafile$str1))
                 rafile$str1 = gsub('0', '-', gsub('1', '+', rafile$str1))
             
             if (is.character(rafile$str2) | is.factor(rafile$str2))
                 rafile$str2 = gsub('0', '-', gsub('1', '+', rafile$str2))
             
             if (is.numeric(rafile$str1))
                 rafile$str1 = ifelse(rafile$str1>0, '+', '-')

             if (is.numeric(rafile$str2))
                 rafile$str2 = ifelse(rafile$str2>0, '+', '-')
             
             rafile$rowid = 1:nrow(rafile)

             bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0
             
             rafile = rafile[which(!bad.ix), ]
             
             if (nrow(rafile)==0)
                 return(GRanges())
             
             seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                 data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

             if (chr.convert)
                 seg$chr = gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr)))
             
             out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
             out = split(out, out$ra.index)
         }
     else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
         {
             ra1 = gr.flip(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
             ra2 = gr.flip(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
             out = grl.pivot(GRangesList(ra1, ra2))             
         }
     
     
     if (keep.features)
         values(out) = rafile[, ]

     return(out)
 }

grl.pivot <-
  function(x)
  {
    if (length(x) == 0)
      return(GRangesList(GRanges(seqlengths = seqlengths(x)), GRanges(seqlengths = seqlengths(x))))
    return(split(unlist(x), rep(1:length(x[[1]]), length(x))))
  }


read_hg <-
function(hg19 = T, fft = F)
  {
    if (fft)
      return(readRDS(REFGENE.FILE.HG19.FFT))
    else
      {
        require(BSgenome)
        if (hg19)
          library(BSgenome.Hsapiens.UCSC.hg19)
        else
          library(BSgenome.Hsapiens.UCSC.hg18)
      }
    return(Hsapiens)
  }

hg_seqlengths <-
function(hg19 = T, chr = F, include.junk = F)
  {
    require(BSgenome)
    hg = read_hg(hg19)

    sl = seqlengths(hg)

    if (!include.junk)
      sl = sl[c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY', 'chrM')]
    
    if (!chr)
      names(sl) = gsub('chr', '', names(sl))

    return(sl)          
  }

grl.unlist <-
function(grl)
  {
    if (length(grl) == 0) ## JEREMIAH
      return(GRanges())
#      return(grl) 
    names(grl) = NULL

    as.df = as.data.frame(grl)
    
    el = as.df$element
    if (is.null(el))
        el = as.df$group
       
    out = unlist(grl)
    out$grl.ix = el
    tmp = rle(el)
    out$grl.iix = unlist(sapply(tmp$lengths, function(x) 1:x))
    values(out) = cbind(values(grl)[out$grl.ix, , drop = FALSE], values(out))
    return(out)
  }

grdt <-
function(x)
 {
      require(data.table)

      ## new approach just directly instantiating data table
      cmd = 'data.table(';
      if (is(x, 'GRanges'))
          {
              was.gr = TRUE
              f = c('seqnames', 'start', 'end', 'strand', 'width')
              f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand', 'as.numeric(width')
              cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
              value.f = names(values(x))              
          }
      else          
          {
              was.gr = FALSE
              value.f = names(x)
          }
      
      if (length(value.f)>0)
          {
              if (was.gr)
                  cmd = paste(cmd, ',', sep = '')
              class.f = sapply(value.f, function(f) eval(parse(text=sprintf("class(x$'%s')", f))))

              .StringSetListAsList = function(x) ### why do I need to do this, bioconductor peeps??
                  {
                      tmp1 = as.character(unlist(x))
                      tmp2 = rep(1:length(x), elementLengths(x))
                      return(split(tmp1, tmp2))                          
                  }

              ## take care of annoying S4 / DataFrame / data.frame (wish-they-were-non-)issues
              as.statement = ifelse(grepl('Integer', class.f), 'as.integer',
                  ifelse(grepl('Character', class.f), 'as.character',
                         ifelse(grepl('StringSetList', class.f), '.StringSetListAsList',
                                ifelse(grepl('StringSet$', class.f), 'as.character',
                                       ifelse(grepl('List', class.f), 'as.list',
                                              ifelse(grepl('List', class.f), 'as.list', 'c'))))))
              cmd = paste(cmd, paste(value.f, '=', as.statement, "(x$'", value.f, "')", sep = '', collapse = ','), sep = '')
          }

      cmd = paste(cmd, ')', sep = '')
      return(eval(parse(text =cmd)))
  }

gr.findoverlaps <-
function(query, subject, ignore.strand = T, first = F,
    qcol = NULL, ## any query meta data columns to add to result
    scol = NULL, ## any subject meta data columns to add to resultx
    max.chunk = 1e13,
    foverlaps = ifelse(is.na(as.logical(Sys.getenv('GRFO_FOVERLAPS'))), TRUE, as.logical(Sys.getenv('GRFO_FOVERLAPS'))) & exists('foverlaps'),
    pintersect = NA,
    verbose = F,
    type = 'any', 
    by = NULL, 
    mc.cores = 1,
    return.type = 'same',
    ...)
  {

     if (type != 'any')
         {
             foverlaps = FALSE
             pintersect = FALSE
         }
      
  if (nchar(foverlaps)==0)
      foverlaps = TRUE

  if (is.na(foverlaps))
      foverlaps = TRUE
  
  isdt <- any(class(query) == 'data.table' )
  if (return.type == 'same')
    return.type <- ifelse(isdt, 'data.table', 'GRanges')
  
  if (!((inherits(subject, 'GRanges') | inherits(subject, 'data.table')) & (inherits(query, 'GRanges') | inherits(query, 'data.table'))))
      stop('both subject and query have to be GRanges or data.table')
  
  if (is.na(pintersect))
    if (isdt)
      pintersect <- length(unique(query$seqnames)) > 50 & length(unique(subject$seqnames)) > 50
    else
      pintersect <- seqlevels(query) > 50 && seqlevels(subject) > 50
  if (is.na(pintersect))
    pintersect <- FALSE

  
  if (!is.null(qcol))
      if (!all(qcol %in% names(values(query))))
          stop('Some qcol are not present in meta data of query')

  if (!is.null(scol))
      if (!all(scol %in% names(values(subject))))
          stop('Some scol are not present in meta data of subject')
          
  if (!is.null(by))
    if (!(by %in% names(values(query)) & by %in% names(values(subject))))
      stop('"by" field must be meta data column of both query and subject')
    
    if ((as.numeric(length(query)) * as.numeric(length(subject))) > max.chunk)
      {
        if (verbose) 
          cat('Overflow .. computing overlaps in chunks.  Adjust max.chunk parameter to gr.findoverlaps to avoid chunked computation\n')
        chunk.size = floor(sqrt(max.chunk));
        ix1 = c(seq(1, length(query), chunk.size), length(query)+1)
        ix2 = c(seq(1, length(subject), chunk.size), length(subject)+1)
        ij = cbind(rep(1:(length(ix1)-1), length(ix2)-1), rep(1:(length(ix2)-1), each = length(ix1)-1))
        if (verbose)
          print(paste('Number of chunks:', nrow(ij)))

        out = do.call('c', mclapply(1:nrow(ij),
            function(x)
                        {
                          if (verbose)
                            cat(sprintf('chunk i = %s-%s (%s), j = %s-%s (%s)\n', ix1[ij[x,1]], ix1[ij[x,1]+1]-1, length(query),
                                        ix2[ij[x,2]], (ix2[ij[x,2]+1]-1), length(subject)))
                          i.chunk = ix1[ij[x,1]]:(ix1[ij[x,1]+1]-1)
                          j.chunk = ix2[ij[x,2]]:(ix2[ij[x,2]+1]-1)
                          out = gr.findoverlaps(query[i.chunk], subject[j.chunk],  ignore.strand = ignore.strand, first = first, pintersect=pintersect, by = by, qcol = qcol, verbose = verbose, foverlaps = foverlaps, scol = scol, type = type, ...)
                          out$query.id = i.chunk[out$query.id]
                          out$subject.id = j.chunk[out$subject.id]
                          return(out)
                        }, mc.cores=mc.cores))

        convert = FALSE
        if ((return.type == 'same' & is(query, 'data.table')) | return.type == 'data.table')
            out = grdt(out)
        return(out)            
      }

  if (foverlaps)
      {
          if (verbose)
              print('overlaps by data.table::foverlaps')
          if (ignore.strand)
              by = c(by, 'seqnames',  'start', 'end')
          else
              by = c(by, 'seqnames', 'strand', 'start', 'end')

          if (!is.data.table(query))
              {
                  names(query) = NULL
                  querydt = grdt(query[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
              }
          else
              {
                  if (!all(by %in% names(query)))
                      stop(paste('the following columns are missing from query:',
                                 paste(by, collapse = ',')))
                      
                  querydt = query[, by, with = FALSE]
              }
          
          if (!is.data.table(subject))
              {
                  names(subject) = NULL
                  subjectdt = grdt(subject[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
              }
          else
              {
                  if (!all(by %in% names(subject)))
                      stop(paste('the following columns are missing from subejct:',
                                 paste(by, collapse = ',')))
                  subjectdt = subject[, by, with = FALSE]
              }
          
          
          ix1 = querydt$query.id = 1:nrow(querydt)
          ix2 = subjectdt$subject.id = 1:nrow(subjectdt)

          querydt = querydt[start<=end, ]
          subjectdt = subjectdt[start<=end, ]
          
          querydt = querydt[, c('query.id', by), with = F]
          subjectdt = subjectdt[, c('subject.id', by), with = F]
          setkeyv(querydt, by)
          setkeyv(subjectdt, by)

         
          h.df = foverlaps(querydt, subjectdt, by.x = by, by.y = by, mult = 'all', type = 'any', verbose = verbose)
          h.df = h.df[!is.na(subject.id) & !is.na(query.id), ]
          h.df[, start := pmax(start, i.start)]
          h.df[, end := pmin(end, i.end)]

          if (verbose)
              cat(sprintf('Generated %s overlaps\n', nrow(h.df)))          
      }
  else
      {

          if (isdt) {
              sn1 <- query$seqnames
              sn2 <- subject$seqnames
          } else {
              sn1 = as.character(seqnames(query))
              sn2 = as.character(seqnames(subject))
          }
          if (is.null(by))
              {
                  ix1 = which(sn1 %in% sn2)
                  ix2 = which(sn2 %in% sn1)
              }
          else
              {
                  by1 = values(query)[, by]
                  by2 = values(subject)[, by]
                  ix1 = which(sn1 %in% sn2 & by1 %in% by2)
                  ix2 = which(sn2 %in% sn1 & by2 %in% by1)
                  by1 = by1[ix1]
                  by2 = by2[ix2]
              }
          
          query.ix = query[ix1]
          subject.ix = subject[ix2]
          sn1 = sn1[ix1]
          sn2 = sn2[ix2]
          
          
          if (pintersect)
              {
                  if (verbose)
                      print('overlaps by pintersect')
                  require(data.table)
                  if (length(sn1)>0 & length(sn2)>0)
                      {

                          if (is.null(by))
                              {
                                  dt1 <- data.table(i=seq_along(sn1), sn=sn1, key="sn")
                                  dt2 <- data.table(j=seq_along(sn2), sn=sn2, key="sn")                
                                  ij <- merge(dt1, dt2, by = 'sn', allow.cartesian=TRUE)
                              }
                          else
                              {
                                  dt1 <- data.table(i=seq_along(sn1), sn=sn1, by = by1, key=c("sn", "by"))
                                  dt2 <- data.table(j=seq_along(sn2), sn=sn2, by = by2, key=c("sn", "by"))
                                  ij <- merge(dt1, dt2, by = c('sn', 'by'), allow.cartesian=TRUE)
                              }

                          if (ignore.strand && isdt)
                              subject$strand <- '*'
                          else if (ignore.strand)
                              strand(subject) = '*'

                          qr <- query.ix[ij$i]
                          sb <- subject.ix[ij$j]
                          if (!isdt) {
                              seqlengths(qr) <- rep(NA, length(seqlengths(qr)))
                              seqlengths(sb) <- rep(NA, length(seqlengths(sb)))
                          }
                          
                          if (!isdt && any(as.character(seqnames(qr)) != as.character(seqnames(sb))))
                              warning('gr.findoverlaps: violated pintersect assumption')

                          ## changed to ranges(qr) etc rather than just GRanges call. Major problem if too many seqlevels
                          if (isdt) {
                              rqr <- IRanges(start=qr$start, end=qr$end)
                              rsb <- IRanges(start=sb$start, end=sb$end)              
                          } else {
                              rqr <- ranges(qr)
                              rsb <- ranges(sb)              
                          }
                          tmp <- pintersect(rqr, rsb, resolve.empty = 'start.x', ...)
                          names(tmp) = NULL
                          non.empty = which(width(tmp)!=0)
                          h.df = as.data.frame(tmp[non.empty])
                          if (isdt)
                              h.df$seqnames <- qr$seqnames[non.empty]
                          else
                              h.df$seqnames <- as.character(seqnames(qr))[non.empty]
                          h.df$query.id = ij$i[non.empty]
                          h.df$subject.id = ij$j[non.empty]
                      }
                  else
                      h.df = data.frame()
              }
          else
              {
                  if (verbose)
                      print('overlaps by findOverlaps')
                  if (isdt) {
                      rqr <- IRanges(start=query.ix$start, end=query.ix$end)
                      rsb <- IRanges(start=subject.ix$start, end=subject.ix$end)              
                  } else {
                      rqr <- ranges(query.ix)
                      rsb <- ranges(subject.ix)              
                  }

                  h <- findOverlaps(rqr, rsb, type = type)
                  r <- ranges(h, rqr, rsb)   
                  h.df <- data.frame(start = start(r), end = end(r), query.id = queryHits(h), subject.id = subjectHits(h), stringsAsFactors = F);
                                        #        sn.query = as.character(seqnames(query))[h.df$query.id]        
                                        #        sn.subject = as.character(seqnames(subject))[h.df$subject.id]
                  sn.query <- sn1[h.df$query.id]
                  sn.subject <- sn2[h.df$subject.id]

                  if (is.null(by))
                      keep.ix <- sn.query == sn.subject
                  else
                      {
                          by.query <- by1[h.df$query.id]
                          by.subject <- by2[h.df$subject.id]
                          keep.ix <- sn.query == sn.subject & by.query == by.subject
                      }
                  
                  h.df <- h.df[keep.ix, ]
                  h.df$seqnames <- sn.query[keep.ix];
              }

          if (!ignore.strand)
              {
                  h.df$strand <- str.query <- as.character(strand(query)[ix1[h.df$query.id]])
                  str.subject <- as.character(strand(subject)[ix2[h.df$subject.id]])
                  h.df <- h.df[which(str.query == str.subject | str.query == '*' | str.subject == '*'),]            
              }
          else if (nrow(h.df)>0)
              h.df$strand = '*'
      }
    
    if (first)
      h.df = h.df[!duplicated(h.df$query.id), ]

     if (return.type=='GRanges') 
       if (nrow(h.df)>0)           
           {
               if (('strand' %in% names(h.df)))
                   out.gr = GRanges(h.df$seqnames, IRanges(h.df$start, h.df$end),
                       query.id = ix1[h.df$query.id], subject.id = ix2[h.df$subject.id], strand = h.df$strand, seqlengths = seqlengths(query))
               else
                   out.gr = GRanges(h.df$seqnames, IRanges(h.df$start, h.df$end),
                       query.id = ix1[h.df$query.id], subject.id = ix2[h.df$subject.id], seqlengths = seqlengths(query))
                   
               if (!is.null(qcol))
                   values(out.gr) = cbind(values(out.gr), values(query)[out.gr$query.id, qcol, drop = FALSE])

               if (!is.null(scol))
                   values(out.gr) = cbind(values(out.gr), values(subject)[out.gr$subject.id, scol, drop = FALSE])               

               return(out.gr)
           }
       else
         return(GRanges(seqlengths = seqlengths(query)))
     else 
         if (nrow(h.df)>0) {
             
             if (!is.data.table(h.df))
                 h.df = as.data.table(h.df)
             h.df$query.id <- ix1[h.df$query.id]
             h.df$subject.id <- ix2[h.df$subject.id]

             if (!is.null(qcol))
                 h.df = cbind(h.df, as.data.table(as.data.frame(values(query))[h.df$query.id, qcol, drop = FALSE]))
             
             if (!is.null(scol))
                 h.df = cbind(h.df, as.data.table(as.data.frame(values(subject))[h.df$subject.id, scol, drop = FALSE]))
             
             if ('i.start' %in% colnames(h.df))
                 h.df[, i.start := NULL]
                          
             if ('i.end' %in% colnames(h.df))
                 h.df[, i.end := NULL]
             
             return(h.df)
       } else {
         return(data.table())
       }
   }


library(optparse)

option_list = list(
    make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input SV VCF file"),
    make_option(c("-o", "--output"), type = "character", default = "circos",  help = "Output basename of pdf to write the graph"),
    make_option(c("-g", "--genes"), type = "logical", default = TRUE,  help = "Add genes to the plot?"),
    make_option(c("-H", "--height"), type = "numeric", default = 10,  help = "Height"),
    make_option(c("-W", "--width"), type = "numeric", default = 10,  help = "Width")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))

bks <- ra_breaks(opt$input)

## filter out discorant only
bks <- bks[-which(mcols(bks)$EVDNC == "DSCRD")]

require(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude <- NULL;
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);

b <- grl.unlist(bks)

## get the gene label dat
basedir <- '/xchip/gistic/Jeremiah/tracks'  
genes <- readRDS(file.path(basedir, 'gr.allgenes.rds'))
genes <- genes[width(genes) < 2e6]
fo <- gr.findoverlaps(b+10e3, genes)
fo <- fo[!duplicated(fo$subject.id)]

gene.dat = data.frame(Chromosome = seqnames(genes[fo$subject.id]), chromStart=start(genes[fo$subject.id]),
  chromEnd=end(genes[fo$subject.id]), Gene=genes$gene[fo$subject.id])
print(gene.dat)

gename.col <- 4;
side <- "in";
track.num <- 1;


x1 = b$grl.iix==1
x2 = b$grl.iix==2
links = data.frame(Chromosome=seqnames(b[x1]), chromStart=start(b[x1]), chromEnd=end(b[x1]), Chromosome.1=seqnames(b[x2]), chromStart.1=start(b[x2]), chromeEnd.1=end(b[x2]))

## plot the PDF
pdf(file=paste0(opt$output,".pdf"), height=opt$height, width=opt$width, compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();
if (opt$genes && nrow(gene.dat) > 0) {
  RCircos.Gene.Connector.Plot(gene.dat, track.num, side);
  track.num <- 2;
  name.col <- 4;
  RCircos.Gene.Name.Plot(gene.dat, name.col,track.num, side);
}
if (nrow(links) > 0)
  RCircos.Link.Plot(links, track.num, by.chromosome=TRUE) ## by.chromosome is for color
dev.off()
