#!/bin/bash
meas=400
databasedir="/cl_tmp/07schadh/su3runs/su3pp/bin"
fraction=0.85
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do 
	echo "#\!/bin/bash" > c40x${t}.sh
        echo "#" >> c40x${t}.sh
        echo "#$ -cwd" >> c40x${t}.sh
        echo "#$ -S /bin/csh" >> c40x${t}.sh
        echo "#$ -V" >> c40x${t}.sh
        echo "#$ -q all.q" >> c40x${t}.sh
        echo "#$ -N analyze${t}" >> c40x${t}.sh
	echo "cd ${PWD}/40x${t}" >> c40x${t}.sh
	echo "/usr/bin/time -o ../40x${t}.log --append ../findcluster.x -s 40 -t $t -n $meas -f $fraction pollev_40x${t}_b6.20.list -b > ../40x${t}.log" >> c40x${t}.sh
	echo "exit 0" >> c40x${t}.sh
done

for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do
	if [ ! -d 40x${t} ]
		then
	        mkdir 40x${t}
	fi
	cd 40x${t}
	pwd
	rm -f pollev_40x${t}_b6.20.list
	ls ${databasedir}/40x${t}/pollev_*m* > pollev_40x${t}_b6.20.list
	cd ..
done
