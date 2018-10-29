#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
sed -n '/R'{'length(R)+1'}' = \[/p' test"$j$i".m > Rs"$j$i".Rdata
echo $j$i' done'
done
done
