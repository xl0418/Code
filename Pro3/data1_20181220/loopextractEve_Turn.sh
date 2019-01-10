#!/bin/bash
for j in {1..5};
do
for i in {1..5};
do 
sed -n '/events = \[\r/,/\];/p' Lphi"$j"psi"$i".m > Levent"$j$i".Rdata
sed -n '/turnover/p' Lphi"$j"psi"$i".m > Lturnover"$j$i".Rdata
A=$(cat Levent"$j$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' Levent"$j$i".Rdata
sed -i 's/%/#/' Lturnover"$j$i".Rdata
echo L$j$i' done'


sed -n '/events = \[\r/,/\];/p' Mphi"$j"psi"$i".m > Mevent"$j$i".Rdata
sed -n '/turnover/p' Mphi"$j"psi"$i".m > Mturnover"$j$i".Rdata
A=$(cat Mevent"$j$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' Mevent"$j$i".Rdata
sed -i 's/%/#/' Mturnover"$j$i".Rdata
echo L$j$i' done'


sed -n '/events = \[\r/,/\];/p' Hphi"$j"psi"$i".m > Hevent"$j$i".Rdata
sed -n '/turnover/p' Hphi"$j"psi"$i".m > Hturnover"$j$i".Rdata
A=$(cat Hevent"$j$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's/\[/structure(c(/' -e 's/\];/),.Dim=c(6,'$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' Hevent"$j$i".Rdata
sed -i 's/%/#/' Hturnover"$j$i".Rdata
echo L$j$i' done'
done
done
