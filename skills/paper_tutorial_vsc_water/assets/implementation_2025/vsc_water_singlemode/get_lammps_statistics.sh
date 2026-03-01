#!/bin/bash

for dir in size_*/
do
	echo "working under $dir"
	cd $dir

	grep -i "Performance:" driver.o* | awk '{print $6}' > lammps_count.txt
	
	cd ..
done	
