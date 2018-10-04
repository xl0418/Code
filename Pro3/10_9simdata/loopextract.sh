#!/bin/bash
for i in {0..1};
do 
sed -n '/events = \[\r/,/\];/p' test0"$i".m > events"$i".Rdata
echo '1 done'
sed -n '/turnover/p' test0"$i".m > turnover"$i".Rdata
echo '2 done'
A=$(cat events"$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' events"$i".Rdata
echo '3 done'
sed -i 's/%/#/' turnover"$i".Rdata
done
