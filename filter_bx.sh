#!/bin/bash

gawk -v cutoff="$1" 'BEGIN {FS = "\t"} { 
     if (/^#/) {
      print $0
     } else {
        match($8, /.*?BX=(.*?);E/, ary);
        split(ary[1], bx, ",")
        if (length(bx) >= cutoff)
          print $0
     }
   }' $2
