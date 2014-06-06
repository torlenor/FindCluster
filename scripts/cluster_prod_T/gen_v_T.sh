#!/bin/bash
for l in `ls ../../pconf/` ; do 
	echo "Processing $l..."
	mkdir $l
	cd $l
	ls ../../../pconf/$l/32x8/pconf.lat.* > pollinput.list
	ln -s ../start_v_T.sh .
	cd ..
done
