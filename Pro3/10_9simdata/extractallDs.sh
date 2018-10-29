#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
sed -n '/D'{'length(D)+1'}' = \[\r/,/\];/p' test"$j$i".m > Ds"$j$i".Rdata
echo $j$i' done'
done
done
