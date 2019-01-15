#!/bin/bash

sed -n '/D'{'length(D)+1'}'/,/\];/p' neutral1e9.m > neuDs.Rdata
echo 'neuDs done'

