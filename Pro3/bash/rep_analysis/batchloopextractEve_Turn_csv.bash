#!/bin/bash
ticks='1e+07'

for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
sed -n '/events = \[\r/,/\];/p' "$ticks"/spatialpara"$ticks"L"$j$i"/Lpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"L"$j$i"/Levent"$j$i"rep"$rep".csv
sed -n '/turnover/p' "$ticks"/spatialpara"$ticks"L"$j$i"/Lpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"L"$j$i"/Lturnover"$j$i"rep"$rep".Rdata
sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' "$ticks"/spatialpara"$ticks"L"$j$i"/Levent"$j$i"rep"$rep".csv
sed -i 's/%/#/' "$ticks"/spatialpara"$ticks"L"$j$i"/Lturnover"$j$i"rep"$rep".Rdata
echo 'Extract turnover and events table to Levent'$j$i'rep'$rep'.csv and Lturnover'$j$i'rep'$rep'.Rdata done'


sed -n '/events = \[\r/,/\];/p' "$ticks"/spatialpara"$ticks"M"$j$i"/Mpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"M"$j$i"/Mevent"$j$i"rep"$rep".csv
sed -n '/turnover/p' "$ticks"/spatialpara"$ticks"M"$j$i"/Mpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"M"$j$i"/Mturnover"$j$i"rep"$rep".Rdata

sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' "$ticks"/spatialpara"$ticks"M"$j$i"/Mevent"$j$i"rep"$rep".csv
sed -i 's/%/#/' "$ticks"/spatialpara"$ticks"M"$j$i"/Mturnover"$j$i"rep"$rep".Rdata
echo 'Extract turnover and events table to Mevent'$j$i'rep'$rep'.csv and Mturnover'$j$i'rep'$rep'.Rdata done'


sed -n '/events = \[\r/,/\];/p' "$ticks"/spatialpara"$ticks"H"$j$i"/Hpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"H"$j$i"/Hevent"$j$i"rep"$rep".csv
sed -n '/turnover/p' "$ticks"/spatialpara"$ticks"H"$j$i"/Hpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"H"$j$i"/Hturnover"$j$i"rep"$rep".Rdata
sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' "$ticks"/spatialpara"$ticks"H"$j$i"/Hevent"$j$i"rep"$rep".csv
sed -i 's/%/#/' "$ticks"/spatialpara"$ticks"H"$j$i"/Hturnover"$j$i"rep"$rep".Rdata
echo 'Extract turnover and events table to Hevent'$j$i'rep'$rep'.csv and Hturnover'$j$i'rep'$rep'.Rdata done'
done
done
done