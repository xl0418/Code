#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
echo $j$i
unix2dos Lphi"$j"psi"$i".m
unix2dos Mphi"$j"psi"$i".m
unix2dos Hphi"$j"psi"$i".m
done
done
