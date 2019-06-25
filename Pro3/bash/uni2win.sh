#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
echo $j$i
unix2dos Lpsi"$j"s_phi"$i".m
unix2dos Mpsi"$j"s_phi"$i".m
unix2dos Hpsi"$j"s_phi"$i".m
done
done
