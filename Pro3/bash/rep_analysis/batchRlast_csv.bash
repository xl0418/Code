#!/bin/bash
ticks='1e+07'
for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
A=$(grep -c 'R'{'length(R)+1'}' = \[' "$ticks"/spatialpara"$ticks"L"$j$i"/LRs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' "$ticks"/spatialpara"$ticks"L"$j$i"/LRs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"L"$j$i"/LRt"$j$i"rep"$rep".Rdata
sed -n '/R = \[/,/\];/p' "$ticks"/spatialpara"$ticks"L"$j$i"/LRt"$j$i"rep"$rep".Rdata > "$ticks"/spatialpara"$ticks"L"$j$i"/LR"$j$i"rep"$rep".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"L"$j$i"/LR"$j$i"rep"$rep".csv
echo 'Extract the last R to LR'$j$i'rep'$rep'.csv done'


A=$(grep -c 'R'{'length(R)+1'}' = \[' "$ticks"/spatialpara"$ticks"M"$j$i"/MRs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' "$ticks"/spatialpara"$ticks"M"$j$i"/MRs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"M"$j$i"/MRt"$j$i"rep"$rep".Rdata
sed -n '/R = \[/,/\];/p' "$ticks"/spatialpara"$ticks"M"$j$i"/MRt"$j$i"rep"$rep".Rdata > "$ticks"/spatialpara"$ticks"M"$j$i"/MR"$j$i"rep"$rep".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"M"$j$i"/MR"$j$i"rep"$rep".csv
echo 'Extract the last R to MR'$j$i'rep'$rep'.csv done'


A=$(grep -c 'R'{'length(R)+1'}' = \[' "$ticks"/spatialpara"$ticks"H"$j$i"/HRs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/R'{'length(R)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/R'{'length(R)+1'}'/R/;:a;n;ba' "$ticks"/spatialpara"$ticks"H"$j$i"/HRs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"H"$j$i"/HRt"$j$i"rep"$rep".Rdata
sed -n '/R = \[/,/\];/p' "$ticks"/spatialpara"$ticks"H"$j$i"/HRt"$j$i"rep"$rep".Rdata > "$ticks"/spatialpara"$ticks"H"$j$i"/HR"$j$i"rep"$rep".csv
sed -i -e 's/R = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"H"$j$i"/HR"$j$i"rep"$rep".csv

echo 'Extract the last R to HR'$j$i'rep'$rep'.csv done'
done
done
done