#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
A=$(grep -c 'D'{'length(D)+1'}' = \[' LDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' LDs"$j$i".Rdata>LDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' LDt"$j$i".Rdata > LD"$j$i".Rdata
C=$(cat LD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/D = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $C)','$(echo $C)'))/' -e 's/;/ /g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $D)'{s/^/,/g}' LD"$j$i".Rdata
sed -i 's/,)/)/g' LD"$j$i".Rdata
echo L$j$i' done'



A=$(grep -c 'D'{'length(D)+1'}' = \[' MDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' MDs"$j$i".Rdata>MDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' MDt"$j$i".Rdata > MD"$j$i".Rdata
C=$(cat MD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/D = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $C)','$(echo $C)'))/' -e 's/;/ /g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $D)'{s/^/,/g}' MD"$j$i".Rdata
sed -i 's/,)/)/g' MD"$j$i".Rdata
echo M$j$i' done'



A=$(grep -c 'D'{'length(D)+1'}' = \[' HDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' HDs"$j$i".Rdata>HDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' HDt"$j$i".Rdata > HD"$j$i".Rdata
C=$(cat HD"$j$i".Rdata |wc -l)
C=$[$C-2]
D=$[$C+1]

sed -i -e 's/D = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $C)','$(echo $C)'))/' -e 's/;/ /g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $D)'{s/^/,/g}' HD"$j$i".Rdata
sed -i 's/,)/)/g' HD"$j$i".Rdata
echo H$j$i' done'
done
done
