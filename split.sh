#!/bin/sh

file=$1
n=`cat $file | wc -l`
line_num=`seq 1 500 $n`
csplit $file `echo $line_num` -f gene
