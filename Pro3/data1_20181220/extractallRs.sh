#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
sed -n '/R'{'length(R)+1'}'/,/\];/p' Lphi"$j"psi"$i".m > LRs"$j$i".Rdata
sed -i 's/ ;//g' LRs"$j$i".Rdata 
echo "Lphi"$j"psi"$i' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' Mphi"$j"psi"$i".m > MRs"$j$i".Rdata
sed -i 's/ ;//g' MRs"$j$i".Rdata 
echo "Mphi"$j"psi"$i' done'
sed -n '/R'{'length(R)+1'}'/,/\];/p' Hphi"$j"psi"$i".m > HRs"$j$i".Rdata
sed -i 's/ ;//g' HRs"$j$i".Rdata 
echo "Hphi"$j"psi"$i' done'
done
done
