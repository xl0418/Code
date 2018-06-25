#!/bin/bash
for i in {1..4};
do 
sed -n '/events = \[\r/,/\];/p' Uloop"$i".m > events"$i".Rdata
sed -n '/turnover/p' Uloop"$i".m > turnover"$i".Rdata
A=$(cat events"$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' events"$i".Rdata
sed -i 's/%/#/' turnover"$i".Rdata
done
