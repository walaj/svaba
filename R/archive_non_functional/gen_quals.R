require(data.table)

f <- fread("grep -v ^# /broad/broadsv/NA12878/GCAT/gcat_illumina_150x/snowman/v115/v115.snowman.indel.vcf", sep='\t')
hdr <- fread("grep ^# /broad/broadsv/NA12878/GCAT/gcat_illumina_150x/snowman/v115/v115.snowman.indel.vcf", sep="\t", header=FALSE)

## rescale the quality score
f[, SL := gsub(".*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:.*?:(.*?)","\\1", V10)]
f[, V6 := SL]

p <- c(hdr$V1, f[, paste(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,sep="\t")])
writeLines(p, "/broad/hptmp/jwala/test.vcf")
