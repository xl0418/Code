#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
sed -n '/events = \[\r/,/\];/p' Lpsi"$j"s_phi"$i".m > Levent"$j$i".csv
sed -n '/turnover/p' Lpsi"$j"s_phi"$i".m > Lturnover"$j$i".Rdata
sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' Levent"$j$i".csv
sed -i 's/%/#/' Lturnover"$j$i".Rdata
echo 'Extract turnover and events table to Levent'$j$i'.csv and Lturnover'$j$i'.Rdata done'


sed -n '/events = \[\r/,/\];/p' Mpsi"$j"s_phi"$i".m > Mevent"$j$i".csv
sed -n '/turnover/p' Mpsi"$j"s_phi"$i".m > Mturnover"$j$i".Rdata

sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' Mevent"$j$i".csv
sed -i 's/%/#/' Mturnover"$j$i".Rdata
echo 'Extract turnover and events table to Mevent'$j$i'.csv and Mturnover'$j$i'.Rdata done'


sed -n '/events = \[\r/,/\];/p' Hpsi"$j"s_phi"$i".m > Hevent"$j$i".csv
sed -n '/turnover/p' Hpsi"$j"s_phi"$i".m > Hturnover"$j$i".Rdata
sed -i -e 's/events = \[//' -e 's/\];//' -e 's/  //g' -e 's/ /,/g' -e 's/;//g' Hevent"$j$i".csv
sed -i 's/%/#/' Hturnover"$j$i".Rdata
echo 'Extract turnover and events table to Hevent'$j$i'.csv and Hturnover'$j$i'.Rdata done'
done
done
