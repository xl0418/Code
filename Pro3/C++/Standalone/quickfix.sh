#!/bin/bash
if [ $# -eq 0 ] ; then
	echo "Please supply the path to the data files"
	exit 1	
fi
if [ -f "$1" ] ; then
	mv "$1" "$1.tmp"
	sed 's/^cpp_skip_label_read_back/\% cpp_skip_label_read_back/g' "$1.tmp" > "$1"
	rm -f "$1.tmp"
else
	cd "$1"      
	for file in *.m; do
		mv $file $file.tmp
		sed 's/^cpp_skip_label_read_back/\% cpp_skip_label_read_back/g' $file.tmp > $file
		rm -f $file.tmp
	done
fi

