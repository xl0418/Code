#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
A=$(grep -c 'R'{'length(R)+1'}' = \[' LRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' LRs"$j$i".Rdata>LRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' LRt"$j$i".Rdata > LR"$j$i".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' LR"$j$i".csv
echo 'Extract the last R to LR'$j$i'.csv done'


A=$(grep -c 'R'{'length(R)+1'}' = \[' MRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' MRs"$j$i".Rdata>MRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' MRt"$j$i".Rdata > MR"$j$i".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' MR"$j$i".csv

echo 'Extract the last R to MR'$j$i'.csv done'


A=$(grep -c 'R'{'length(R)+1'}' = \[' HRs"$j$i".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' HRs"$j$i".Rdata>HRt"$j$i".Rdata
sed -n '/R = \[/,/\];/p' HRt"$j$i".Rdata > HR"$j$i".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' HR"$j$i".csv

echo 'Extract the last R to HR'$j$i'.csv done'
done
done
