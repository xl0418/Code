#!/bin/bash
ticks='1e+07'

for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
sed -n '/D'{'length(D)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"L"$j$i"/Lpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"L"$j$i"/LDs"$j$i"rep"$rep".Rdata
echo 'Extract all Ds in Lpsi'$j's_phi'$i'rep'$rep' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"M"$j$i"/Mpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"M"$j$i"/MDs"$j$i"rep"$rep".Rdata
echo 'Extract all Ds in Mpsi'$j's_phi'$i'rep'$rep' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' "$ticks"/spatialpara"$ticks"H"$j$i"/Hpsi"$j"s_phi"$i"rep"$rep".m > "$ticks"/spatialpara"$ticks"H"$j$i"/HDs"$j$i"rep"$rep".Rdata
echo 'Extract all Ds in Hpsi'$j's_phi'$i'rep'$rep' done'
done
done
done