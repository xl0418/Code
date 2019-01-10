#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
sed -n '/D'{'length(D)+1'}' = \[\r/,/\];/p' test"$j$i".m > D"$j$i".Rdata
A=$(cat D"$j$i".Rdata |wc -l)
A=$[$A-2]
echo $A
B=$[$A+1]
sed -i -e 's//D'{'length(D)+1'}' = \[/D = structure(c(/' -e 's/\];/),.Dim=c('$(echo $A)','$(echo $A)'))/' -e 's/;//g' -e '2,${s/ /,/g}' -e '2,${s/,,/ /g}' -e '3,'$(echo $B)'{s/^/,/g}' D"$j$i".Rdata
echo $j$i' done'
done
done
