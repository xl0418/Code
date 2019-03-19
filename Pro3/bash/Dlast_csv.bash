#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
A=$(grep -c 'D'{'length(D)+1'}' = \[' HDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' HDs"$j$i".Rdata>HDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' HDt"$j$i".Rdata>HD"$j$i".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' HD"$j$i".csv
echo H$j$i' done'

A=$(grep -c 'D'{'length(D)+1'}' = \[' MDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' MDs"$j$i".Rdata>MDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' MDt"$j$i".Rdata>MD"$j$i".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' MD"$j$i".csv
echo M$j$i' done'

A=$(grep -c 'D'{'length(D)+1'}' = \[' LDs"$j$i".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' LDs"$j$i".Rdata>LDt"$j$i".Rdata
sed -n '/D = \[/,/\];/p' LDt"$j$i".Rdata>LD"$j$i".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' LD"$j$i".csv
echo L$j$i' done'
done
done
