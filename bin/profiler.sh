#!/usr/bin/env bash
"$@" & # Run the given command line in the background.
pid=$! peak=0
rm mem.log

## make the Rscript
rm memplot.R
echo "require(ggplot2); df = read.table('mem.log', header=FALSE); df\$time = seq(from=10,to=10*nrow(df), by=10);" >> memplot.R;
echo "pdf('memgraph.pdf'); g <- ggplot(df, aes(x=time, y=V1/1024)) + geom_line() + xlab('Time (s)') + ylab('Memory (Mb)') + theme_bw(); print(g); dev.off();" >> memplot.R;
chmod +x memplot.R;

while true; do
  sleep 0.1
  sample="$(ps -o rss= $pid 2> /dev/null)" || break
  let peak='sample > peak ? sample : peak'
  echo "$sample" >> mem.log
done
echo "Peak: $peak" 1>&2

Rscript memplot.R
