#!/use/bin/env/ Rscript

require(ggplot2)
require(data.table)

snow <- ra_breaks("/xchip/gistic/Jeremiah/Projects/SnowmanPaper/Benchmark/snow/chr1.broad-snowman.DATECODE.somatic.sv.vcf")

## read it
dt <- fread("/xchip/gistic/Jeremiah/Projects/SnowmanPaper/150805benchmark.csv") 
dt [, mean_cc := mean(contig_coverage), by=c('kmer_corr', 'coverage', 'error_rate')]
dt [, se_cc := sd(contig_coverage), by=c('kmer_corr', 'coverage', 'error_rate')]
setkey(dt, kmer_corr, coverage, error_rate)
dt <- unique(dt)

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
