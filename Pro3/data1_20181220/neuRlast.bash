#!/bin/bash

A=$(grep -c 'R'{'length(R)+1'}' = \[' neuRs.Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' neuRs.Rdata>neuRt.Rdata
sed -n '/R = \[/,/\];/p' neuRt.Rdata > neuRl.Rdata 
C=$(cat neuDl.Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/R = \[/R=structure(c(/' -e 's/\];/),.Dim=c(1,'$(echo $C)'))/' -e 's/ /,/g' neuRl.Rdata
echo 'neuRl done'
