#!/bin/bash

for dir in size_*/
do
	echo "working under $dir"
	cd $dir

	grep -i "on time step" main.o* | tail -n +2  | awk '{print $6}' > time_count.txt
	
	cd ..
done	
