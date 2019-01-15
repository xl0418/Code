#!/bin/bash

sed -n '/R'{'length(R)+1'}'/,/\];/p' neutral1e9.m > neuRs.Rdata
sed -i 's/ ;//g' neuRs.Rdata 
echo "neuR done'

