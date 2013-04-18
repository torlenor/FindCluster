#!/bin/bash
meas=400
radius=0.90
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do 
	echo "#\!/bin/bash" > c40x${t}.sh
	echo "cd 40x${t}" >> c40x${t}.sh
	echo "/usr/bin/time -o ../40x${t}.log --append ../findcluster.x -s 40 -t $t -n $meas -r $radius pollev_40x${t}_b6.20.list -b > ../40x${t}.log" >> c40x${t}.sh
	echo "cd .." >> c40x${t}.sh
	chmod +x c40x${t}.sh
done

for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do
        if [ ! -d 40x${t} ]
                then
                mkdir 40x${t}
        fi  
	cd 40x${t}
	pwd
	rm pollev_40x${t}_b6.20.list
	ls ../../b6.2/40x${t}/pollev_*m* > pollev_40x${t}_b6.20.list
	cd ..
done

rm *.log

parallel -j 4 ::: \
./c40x2.sh \
./c40x3.sh \
./c40x4.sh \
./c40x5.sh \
./c40x6.sh \
./c40x7.sh \
./c40x8.sh \
./c40x9.sh \
./c40x10.sh \
./c40x11.sh \
./c40x12.sh \
./c40x14.sh \
./c40x16.sh \
./c40x18.sh \
./c40x20.sh

rm percc.res
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do cat 40x${t}/npercc.res ; done > percc.res

rm cutrate.res
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do echo "$t `cat 40x${t}.log | grep Cut`" >> cutrate.res ; done
