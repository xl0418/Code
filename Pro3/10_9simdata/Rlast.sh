#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
A=$(grep -c 'R'{'length(R)+1'}' = \[' Rs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' Rs"$j$i".Rdata>Rt"$j$i".Rdata
sed -n '/R = \[/p' Rt"$j$i".Rdata > R"$j$i".Rdata
C=$(cat D"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/R = \[ /R=structure(c(/' -e 's/ \];/),.Dim=c(1,'$(echo $C)'))/' -e 's/ /,/g' R"$j$i".Rdata
echo $j$i' done'
done
done
