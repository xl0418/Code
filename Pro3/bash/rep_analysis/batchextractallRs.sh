#!/bin/bash
ticks='1e+07'

for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
sed -n '/R'{'length(R)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"L"$j$i"/Lpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"L"$j$i"/LRs"$j$i"rep"$rep".Rdata
sed -i 's/ ;//g' "$ticks"/spatialpara"$ticks"L"$j$i"/LRs"$j$i"rep"$rep".Rdata 
echo "Extract all Rs in Lpsi"$j"s_phi"$i'rep'$rep' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"M"$j$i"/Mpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"M"$j$i"/MRs"$j$i"rep"$rep".Rdata
sed -i 's/ ;//g' "$ticks"/spatialpara"$ticks"M"$j$i"/MRs"$j$i"rep"$rep".Rdata 
echo 'Extract all Rs in Mpsi'$j's_phi'$i'rep'$rep' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"H"$j$i"/Hpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"H"$j$i"/HRs"$j$i"rep"$rep".Rdata
sed -i 's/ ;//g' "$ticks"/spatialpara"$ticks"H"$j$i"/HRs"$j$i"rep"$rep".Rdata 
echo 'Extract all Rs in Hpsi'$j's_phi'$i'rep'$rep' done'
done
done
done