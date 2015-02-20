#!/bin/bash  
for x in $*; do  
		filename=${x##*/}
		header_file="header1.txt"
		orig_header_file="header.txt"
		filestring="/* File: ${filename} */"
		echo "adding to file = $filename"
		echo "$filestring" > $header_file
		cat $orig_header_file >> $header_file
		HEADERLEN=`wc -l ${header_file} | awk '{print $1}' `
		head -n $HEADERLEN $x | diff ${header_file} - || ( ( cat ${header_file}; echo; cat $x) > /tmp/file;  
				mv /tmp/file $x )
done 
