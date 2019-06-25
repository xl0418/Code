#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
sed -n '/D'{'length(D)+1'}'/,/\];/p' Lpsi"$j"s_phi"$i".m > LDs"$j$i".Rdata
echo 'Extract all Ds in Lpsi'$j's_phi'$i' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' Mpsi"$j"s_phi"$i".m > MDs"$j$i".Rdata
echo 'Extract all Ds in Mpsi'$j's_phi'$i' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' Hpsi"$j"s_phi"$i".m > HDs"$j$i".Rdata
echo 'Extract all Ds in Hpsi'$j's_phi'$i' done'
done
done
