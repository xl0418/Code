#!/bin/bash
for j in {1..5};
do
for i in {1..5};
do 
sed -n '/D'{'length(D)+1'}'/,/\];/p' Lphi"$j"psi"$i".m > LDs"$j$i".Rdata
echo 'Lphi'$j'psi'$i' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' Mphi"$j"psi"$i".m > MDs"$j$i".Rdata
echo 'Mphi'$j'psi'$i' done'
sed -n '/D'{'length(D)+1'}'/,/\];/p' Hphi"$j"psi"$i".m > HDs"$j$i".Rdata
echo 'Hphi'$j'psi'$i' done'
done
done
