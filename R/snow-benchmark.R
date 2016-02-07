#!/use/bin/env/ Rscript

require(ggplot2)
require(data.table)

.libPaths = c("/xchip/gistic/Jeremiah/R", "/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/lib64/R/library")

library(optparse)

AVAIL_MODES <- c("realign-test")

option_list = list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input SV VCF file"),
  make_option(c("-o", "--output"), type = "character", default = "no_id",  help = "Output annotation name"),
  make_option(c("-m", "--mode"), type = "character", default = "realign-test",  help = "Benchmarking mode to analyze")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$input))
  stop(print_help(parseobj))


if (!opt$mode %in% AVAIL_MODES)
  stop(paste("Mode", opt$mode, "must be one of", paste(AVAIL_MODES, collapse=", ")))

.format_flag_df <- function(dfh2) {
  dfh2 <- data.table(dfh[dfh$del_size == "0" & dfh$snv_rate == "0" & dfh$ins_size != "0", ])
  dfh2$index = seq(nrow(dfh2))
  dfh2[, total := sum(c(wrong_align, too_align, no_align, correct)), by=index]
  dfh2 = data.frame(dfh2)
  dfh2[, c("no_align", "wrong_align", "too_align","correct")] = sweep(dfh2[, c('no_align', 'wrong_align', 'too_align', 'correct')], 1, dfh2$total, "/")
  dfh2 <- melt(dfh2, id.vars=c("width","ins_size", "del_size", "snv_rate", "total", "index"))
  return (dfh2)
}



if (opt$mode == "realign-test") {

  ff <- fread(opt$input, header=TRUE)

  df <- data.frame();
  dfh <- data.frame();
  for (j in unique(ff$width)) {
    for (k in unique(ff$snv_rate)) {
      for (d in unique(ff$del_size)) {
        for (ii in unique(ff$ins_size)) {
          #ii = ifelse(d == 0, 0, ii) ## if del is zero, then make ins zero too, since just want to evaluate zeros
          #d  = ifelse(ii == 0, 0, d) ## if ins is zero, then make del zero too, since we just want to evaulate zeros
          if (ii == 0 || d == 0) {
            ix  <- ff$width == j & ff$snv_rate == k & ff$ins_size == d & ff$del_size == ii
            if (sum(ix)) {
              cf <- ecdf(ff$num_align[ix])
              s <- seq(0,max(ff$num_align))
              df <- rbind(df, data.frame(num_align=s, cdf=cf(s), width=as.character(j), snv_rate=as.character(k), del_size = as.character(d), ins_size=as.character(ii), one.x = 0, one.y = cf(1)))

              dfh <- rbind(dfh, data.frame(no_align=sum(ix & ff$num_align == 0),
                                           wrong_align=sum(ix & ff$correct_hit_num == -1 & ff$num_align == 1),
                                           too_align=sum(ix & ff$num_align > 1),
                                           correct=sum(ix & ff$correct_hit_num == 0),
                                           width=as.character(j),
                                           del_size = as.character(d),
                                           ins_size = as.character(ii),
                                           snv_rate = as.character(k)))
            }
          }
        }
      }
    }
  }

  ## SNV ONLY
  df2 <- df[df$ins_size == "0" & df$del_size == "0", ]
  g <- ggplot(data=df2) + geom_line(aes(x=num_align, y=cdf, color=width)) + geom_point(aes(x=one.x, y=one.y, color=width)) +  theme_bw() + xlab("Number of alignments") + ylab("CDF") + facet_wrap(~ snv_rate, nrow=1) + scale_y_continuous(limits=c(min(df2$one.y)-0.1,1), breaks=seq(0,1,by=0.05)) + labs(color="Sequence Length")
  pdf("~/public_html/realign_test_snv_cigcheck.pdf", width=7, height=2)
  print(g)
  dev.off()

  ## DEL ONLY
  df2 <- df[df$ins_size == "0" & df$snv_rate == "0", ]
  g <- ggplot(data=df2) + geom_line(aes(x=num_align, y=cdf, color=width)) + geom_point(aes(x=one.x, y=one.y, color=width)) +  theme_bw() + xlab("Number of alignments") + ylab("CDF") + facet_wrap(~ del_size, nrow=1) + scale_y_continuous(limits=c(min(df2$one.y)-0.1,1), breaks=seq(0,1,by=0.2)) + labs(color="Sequence Length")
  pdf("~/public_html/realign_test_del_cigcheck.pdf", width=7, height=2)
  print(g)
  dev.off()

  ## INS ONLY
  df2 <- df[df$del_size == "0" & df$snv_rate == "0" & df$ins_size != "0", ]
  g <- ggplot(data=df2) + geom_line(aes(x=num_align, y=cdf, color=width)) + geom_point(aes(x=one.x, y=one.y, color=width)) +  theme_bw() + xlab("Number of alignments") + ylab("CDF") + facet_wrap(~ ins_size, nrow=1) + scale_y_continuous(limits=c(min(df2$one.y)-0.1,1), breaks=seq(0,1,by=0.2)) + labs(color="Sequence Length")
  pdf("~/public_html/realign_test_ins_cigcheck.pdf", width=7, height=2); print(g); dev.off()
  dfh2 <- .format_flag_df(df2)
  g2 <- ggplot(data=dfh2,aes(x=factor(width),y=value,fill=factor(variable))) + geom_bar(position="stack", stat='identity') + facet_wrap(~ ins_size, nrow=1) + scale_fill_manual(values=c("no_align"="black", "wrong_align"="red", "too_align"="purple","correct"="dark green"), labels=c("Unmapped", "Incorrect alignment", "> 1 alignment", "Accurate"), name="Alignment") + scale_y_continuous(breaks=seq(0,1,by=0.25), labels=c("0", "25", "50", "75", "100"), name="Percentage") + xlab("Sequence Length")
  pdf("~/public_html/realign_ins_flag.pdf", width=7, height=2); print(g2); dev.off()
 
  #g <- ggplot(data=ff) + geom_histogram(aes(x=num_aligns)) + scale_y_log10(limits=c(1, 10)) + theme_bw()

  pdf("~/public_html/realign_test.pdf", width=7, height=7)
  print(g)
  dev.off()
}


if (FALSE) {

snow <- ra_breaks("/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/snow2/chr1.broad-snowman.DATECODE.somatic.sv.vcf")

simd <- fread("/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/connections.tsv")
gr.sim <- with(simd, GRanges(c(V1,V1)+1, IRanges(c(V2,V4), width=1), strand=ifelse(c(V3, V5)=='+', '+', '-'), id=rep(seq(nrow(simd)),each=2)))
grl.sim <- split(gr.sim, gr.sim$id)

ra.overlaps(snow, grl.sim, pad=10, ignore.strand=TRUE)

## read it
dt <- fread("/xchip/gistic/Jeremiah/Projects/SnowmanPaper/150805benchmark.csv") 
dt [, mean_cc := mean(contig_coverage), by=c('kmer_corr', 'coverage', 'error_rate')]
dt [, se_cc := sd(contig_coverage), by=c('kmer_corr', 'coverage', 'error_rate')]
setkey(dt, kmer_corr, coverage, error_rate)
dt <- unique(dt)

dt <- readRDS("/xchip/gistic/Jeremiah/tracks/100map.dt.rds")

## relabel names
en <- c("0"="Error Rate: 0", "0.001"="Error Rate: 1e-3", "0.005"="Error Rate 5e-3",
        "0.01"="Error Rate: 0.01", "0.03"="Error Rate: 0.03", "0.05"="Error Rate: 0.05", "0.1"="Error Rate: 0.1")
.labeller <- function(variable,value){
    return(en[value])
  }

## plot it
df = data.frame(dt)
df$error_rate = as.character(df$error_rate)
df$kmer_corr = factor(df$kmer_corr, levels=c(0,1))
g <- ggplot(df) + geom_line(aes(x=coverage, y=mean_cc, color=kmer_corr), size=1) +
  geom_errorbar(aes(x=coverage,ymin=mean_cc-se_cc, ymax=mean_cc+se_cc, color=kmer_corr), width=1) +
  theme_bw() + xlab("Coverage") + ylab("Percent Re-assembled") +
  scale_x_continuous(breaks=seq(0,40,by=5)) + facet_grid(error_rate ~ ., labeller=.labeller) +
  coord_cartesian(ylim=c(0,1.3)) +
  scale_y_continuous(breaks=seq(0,1,by=0.25))
pdf("~/public_html/plot.pdf", width=6, height=12); print(g); dev.off()


}
