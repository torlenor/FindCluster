#!/bin/bash

su2params="-s 8 -t 4 -b 5.00 -n 20 -k 5 -e 40 -c -m"
su3params="-s 8 -t 4 -b 6.00 -n 20 -k 5 -e 40 -c -m"

function cleanup {
	cd ..
  rm -r build_test
}

echo "Compiling..."
mkdir build_test
cd build_test
cmake ../
make -j 2
if [ -f "./bin/findcluster.x" ]
then
	echo "Compile successful !"
  cp ./bin/findcluster.x .
else
	echo "Compile NOT successful !"
  cleanup
	exit
fi
echo
echo "Performing test run..."
cp ../testfiles/* .

/usr/bin/time ./findcluster.x -s 32 -t 8 -n 5 qcd.list -b -a -w -l -r -f 0.00 | tee qcd_status_f0.00.out
/usr/bin/time ./findcluster.x -s 32 -t 8 -n 5 qcd.list -b -a -w -l -r -f 0.45 | tee qcd_status_f0.45.out
mkdir tmp_qcd
mv *.res tmp_qcd
mv qcd_status_f0.*.out tmp_qcd
/usr/bin/time ./findcluster.x -s 40 -t 9 -n 5 puregauge.list -b -a -w -l -r -f 0.00 -o | tee puregauge_status_f0.00.out
/usr/bin/time ./findcluster.x -s 40 -t 9 -n 5 puregauge.list -b -a -w -l -r -f 0.45 -o | tee puregauge_status_f0.45.out
mkdir tmp_puregauge
mv *.res tmp_puregauge
mv puregauge_status_f0.*.out tmp_puregauge

echo
echo "Cleanup in (cancel with ctrl-c) ..."
for ii in `seq 10 -1 1` ; do 
  echo $ii
  sleep 1
done
cleanup
cd ..
