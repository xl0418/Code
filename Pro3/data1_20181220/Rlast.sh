#!/bin/bash
for j in {1..6};
do
for i in {1..6};
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


A=$(grep -c 'R'{'length(R)+1'}' = \[' MRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' MRs"$j$i".Rdata>MRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' MRt"$j$i".Rdata > MR"$j$i".Rdata
C=$(cat MD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/R = \[/R=structure(c(/' -e 's/\];/),.Dim=c(1,'$(echo $C)'))/' -e 's/ /,/g' MR"$j$i".Rdata
echo 'MR'$j$i' done'


A=$(grep -c 'R'{'length(R)+1'}' = \[' HRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' HRs"$j$i".Rdata>HRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' HRt"$j$i".Rdata > HR"$j$i".Rdata
C=$(cat HD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/R = \[/R=structure(c(/' -e 's/\];/),.Dim=c(1,'$(echo $C)'))/' -e 's/ /,/g' HR"$j$i".Rdata
echo 'HR'$j$i' done'
done
done
