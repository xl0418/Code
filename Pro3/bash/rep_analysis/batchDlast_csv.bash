#!/bin/bash
ticks='1e+07'
for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
A=$(grep -c 'D'{'length(D)+1'}' = \[' "$ticks"/spatialpara"$ticks"H"$j$i"/HDs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' "$ticks"/spatialpara"$ticks"H"$j$i"/HDs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"H"$j$i"/HDt"$j$i"rep"$rep".Rdata
sed -n '/D = \[/,/\];/p' "$ticks"/spatialpara"$ticks"H"$j$i"/HDt"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"H"$j$i"/HD"$j$i"rep"$rep".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"H"$j$i"/HD"$j$i"rep"$rep".csv
echo 'Extract the last D to HD'$j$i'rep'$rep'.csv done'

A=$(grep -c 'D'{'length(D)+1'}' = \[' "$ticks"/spatialpara"$ticks"M"$j$i"/MDs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' "$ticks"/spatialpara"$ticks"M"$j$i"/MDs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"M"$j$i"/MDt"$j$i"rep"$rep".Rdata
sed -n '/D = \[/,/\];/p' "$ticks"/spatialpara"$ticks"M"$j$i"/MDt"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"M"$j$i"/MD"$j$i"rep"$rep".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"M"$j$i"/MD"$j$i"rep"$rep".csv
echo 'Extract the last D to MD'$j$i'rep'$rep'.csv done'

A=$(grep -c 'D'{'length(D)+1'}' = \[' "$ticks"/spatialpara"$ticks"L"$j$i"/LDs"$j$i"rep"$rep".Rdata)
B=$[$A-1]
sed '/D'{'length(D)+1'}'/{G;s/\nX\{'$B'\}//;tend;x;s/^/X/;x;P;d};p;d;:end;s/D'{'length(D)+1'}'/D/;:a;n;ba' "$ticks"/spatialpara"$ticks"L"$j$i"/LDs"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"L"$j$i"/LDt"$j$i"rep"$rep".Rdata
sed -n '/D = \[/,/\];/p' "$ticks"/spatialpara"$ticks"L"$j$i"/LDt"$j$i"rep"$rep".Rdata>"$ticks"/spatialpara"$ticks"L"$j$i"/LD"$j$i"rep"$rep".csv
sed -i -e 's/D = \[//' -e 's/\];//' -e 's/ /,/g' -e 's/,;//g' "$ticks"/spatialpara"$ticks"L"$j$i"/LD"$j$i"rep"$rep".csv
echo 'Extract the last D to LD'$j$i'rep'$rep'.csv done'
done
done
done