#!/bin/bash
databasedir="/cl_tmp/07schadh/su3runs/su3pp/bin"
Ns=40
Nt="2 3 4 5 6 7 8 9 10 11 12 14 16 18 20"
meas=200

fraction=0.00

box=0 # 0/1 to disable/enable box counting

for t in $Nt ; do 
	echo "#\!/bin/bash" > c${Ns}x${t}.sh
        echo "#" >> c${Ns}x${t}.sh
        echo "#$ -cwd" >> c${Ns}x${t}.sh
        echo "#$ -S /bin/csh" >> c${Ns}x${t}.sh
        echo "#$ -V" >> c${Ns}x${t}.sh
        echo "#$ -q all.q" >> c${Ns}x${t}.sh
        echo "#$ -N analyze${t}" >> c${Ns}x${t}.sh
	echo "setenv MV2_ENABLE_AFFINITY 0"  >> c${Ns}x${t}.sh
	echo "cd ${PWD}/${Ns}x${t}" >> c${Ns}x${t}.sh
	if [ $box -eq 1 ] ; then
		echo "/usr/bin/time -o ../${Ns}x${t}.log --append ../findcluster.x -s $Ns -t $t -n $meas -f $fraction pollev_${Ns}x${t}.list -b > ../${Ns}x${t}.log" >> c${Ns}x${t}.sh
	else
		echo "/usr/bin/time -o ../${Ns}x${t}.log --append ../findcluster.x -s $Ns -t $t -n $meas -f $fraction pollev_${Ns}x${t}.list > ../${Ns}x${t}.log" >> c${Ns}x${t}.sh
	fi
	echo "exit 0" >> c${Ns}x${t}.sh
done

for t in $Nt ; do
	if [ ! -d ${Ns}x${t} ]
		then
	        mkdir ${Ns}x${t}
	fi
	cd ${Ns}x${t}
	pwd
	rm -f pollev_${Ns}x${t}.list
	ls ${databasedir}/${Ns}x${t}/pollev_*m* > pollev_${Ns}x${t}.list
	cd ..
done
