#!/bin/bash
for j in {1..6};
do
for i in {1..6};
do 
rm HDs"$j$i".Rdata
rm HDt"$j$i".Rdata
rm HRs"$j$i".Rdata
rm HRt"$j$i".Rdata
echo 'Remove HDs/t'$j$i'.Rdata and HRs/t'$j$i'.Rdata done'

rm MDs"$j$i".Rdata
rm MDt"$j$i".Rdata
rm MRs"$j$i".Rdata
rm MRt"$j$i".Rdata
echo 'Remove MDs/t'$j$i'.Rdata and MRs/t'$j$i'.Rdata done'

rm LDs"$j$i".Rdata
rm LDt"$j$i".Rdata
rm LRs"$j$i".Rdata
rm LRt"$j$i".Rdata
echo 'Remove LDs/t'$j$i'.Rdata and LRs/t'$j$i'.Rdata done'

done
done
