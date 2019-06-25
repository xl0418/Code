#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
sed -n '/R'{'length(R)+1'}'/,/\];/p' Lpsi"$j"s_phi"$i".m > LRs"$j$i".Rdata
sed -i 's/ ;//g' LRs"$j$i".Rdata 
echo "Extract all Rs in Lpsi"$j"s_phi"$i' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' Mpsi"$j"s_phi"$i".m > MRs"$j$i".Rdata
sed -i 's/ ;//g' MRs"$j$i".Rdata 
echo "Extract all Rs in Mpsi"$j"s_phi"$i' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' Hpsi"$j"s_phi"$i".m > HRs"$j$i".Rdata
sed -i 's/ ;//g' HRs"$j$i".Rdata 
echo "Extract all Rs in Hpsi"$j"s_phi"$i' done'
done
done
