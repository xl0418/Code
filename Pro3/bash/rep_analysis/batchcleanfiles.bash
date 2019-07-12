#!/bin/bash
ticks='1e+07'

for j in {1..3};
do
for i in {1..6};
do 
for rep in {61..100};
do
rm "$ticks"/spatialpara"$ticks"H"$j$i"/HDs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"H"$j$i"/HDt"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"H"$j$i"/HRs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"H"$j$i"/HRt"$j$i"rep"$rep".Rdata
echo 'Remove HDs/t'$j$i'rep'$rep'.Rdata and HRs/t'$j$i'rep'$rep'.Rdata done'

rm "$ticks"/spatialpara"$ticks"M"$j$i"/MDs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"M"$j$i"/MDt"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"M"$j$i"/MRs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"M"$j$i"/MRt"$j$i"rep"$rep".Rdata
echo 'Remove MDs/t'$j$i'rep'$rep'.Rdata and MRs/t'$j$i'rep'$rep'.Rdata done'

rm "$ticks"/spatialpara"$ticks"L"$j$i"/LDs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"L"$j$i"/LDt"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"L"$j$i"/LRs"$j$i"rep"$rep".Rdata
rm "$ticks"/spatialpara"$ticks"L"$j$i"/LRt"$j$i"rep"$rep".Rdata
echo 'Remove LDs/t'$j$i'rep'$rep'.Rdata and LRs/t'$j$i'rep'$rep'.Rdata done'

done
done
done