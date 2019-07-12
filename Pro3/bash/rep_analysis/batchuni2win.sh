#!/bin/bash
ticks='1e+07'

for j in {1..6};
do
for i in {1..6};
do 
for rep in {61..100};
do
echo $j$i$rep
unix2dos "$ticks"/spatialpara"$ticks"L"$j$i"/Lpsi"$j"s_phi"$i"rep"$rep".m
unix2dos "$ticks"/spatialpara"$ticks"M"$j$i"/Mpsi"$j"s_phi"$i"rep"$rep".m
unix2dos "$ticks"/spatialpara"$ticks"H"$j$i"/Hpsi"$j"s_phi"$i"rep"$rep".m
done
done
done