#!/bin/bash
for j in {0..4};
do
for i in {0..4};
do 
echo $j$i
unix2dos test"$j$i".m
done
done
