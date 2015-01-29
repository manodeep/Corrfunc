#!/bin/bash  
header_file="header.txt"
HEADERLEN=`wc -l ${header_file} | awk '{print $1}' `
totlen=$((HEADERLEN+1))
for x in $*; do  
filename=${x##*/} 
head -n $totlen $x |tail -n  $HEADERLEN | diff ${header_file} - || ( ( echo "/* File: ${filename} */"; cat ${header_file}; echo; cat $x) > /tmp/file;  
mv /tmp/file $x )  
done 
