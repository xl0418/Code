#!/bin/bash
for j in {1..5};
do
for i in {1..5};
do 
A=$(grep -c 'R'{'length(R)+1'}' = \[' LRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' LRs"$j$i".Rdata>LRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' LRt"$j$i".Rdata > LR"$j$i".Rdata
C=$(cat LD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/R = \[/R=structure(c(/' -e 's/\];/),.Dim=c(1,'$(echo $C)'))/' -e 's/ /,/g' LR"$j$i".Rdata
echo 'LR'$j$i' done'
done
done
