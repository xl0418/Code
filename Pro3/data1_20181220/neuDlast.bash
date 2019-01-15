#!/bin/bash

A=$(grep -c 'D'{'length(D)+1'}' = \[' neuDs.Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' neuDs.Rdata>neuDt.Rdata
sed -n '/D = \[/,/\];/p' neuDt.Rdata > neuDl.Rdata
C=$(cat neuDl.Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/D = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $C)','$(echo $C)'))/' -e 's/;/ /g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $D)'{s/^/,/g}' neuDl.Rdata
sed -i 's/,)/)/g' neuDl.Rdata
echo L$j$i' done'


