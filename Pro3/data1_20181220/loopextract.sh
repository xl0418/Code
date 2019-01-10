#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
sed -n '/events = \[\r/,/\];/p' test"$j$i".m > events"$j$i".Rdata
sed -n '/turnover/p' test"$j$i".m > turnover"$j$i".Rdata
A=$(cat events"$j$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' events"$j$i".Rdata
sed -i 's/%/#/' turnover"$j$i".Rdata
echo $j$i' done'
done
done
