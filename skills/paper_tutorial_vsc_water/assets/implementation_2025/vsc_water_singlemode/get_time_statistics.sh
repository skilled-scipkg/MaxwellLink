#!/bin/bash

for dir in size_*/
do
	echo "working under $dir"
	cd $dir

	grep -i "\[SingleModeCavity\] Completed" main.o* | tail -n +2  | awk '{print $7}' > time_count.txt
	
	cd ..
done	
