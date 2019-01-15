#!/bin/bash

sed -n '/events = \[\r/,/\];/p' neutral1e9.m > neuevent.Rdata
sed -n '/turnover/p' neutral1e9.m > neuturnover.Rdata
A=$(cat neuevent.Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' neuevent.Rdata
sed -i 's/%/#/' neuturnover.Rdata
echo L$j$i' done'



