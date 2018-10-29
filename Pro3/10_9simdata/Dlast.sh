#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
A=$(grep -c 'D'{'length(D)+1'}' = \[' Ds"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' Ds"$j$i".Rdata>Dt"$j$i".Rdata
sed -n '/D = \[\r/,/\];/p' Dt"$j$i".Rdata > D"$j$i".Rdata
C=$(cat D"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/D = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $C)','$(echo $C)'))/' -e 's/;/ /g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $D)'{s/^/,/g}' D"$j$i".Rdata
sed -i 's/,)/)/g' D"$j$i".Rdata
echo $j$i' done'
done
done
