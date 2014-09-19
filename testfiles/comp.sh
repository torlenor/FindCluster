#!/bin/bash
cd tmp_qcd
if [ $? -ne 0 ] ; then echo "ERROR: Exiting..." ; exit ; fi
echo "Testing qcd results (if no output, everything OK!)"
for l in *.res ; do
  diff $l ../comp_qcd_wupper_32x8.bin/$l -q
  if [ $? -ne 0 ] ; then
    echo "File $l differs (or is not there)!"
  fi
done
cd ..

cd tmp_puregauge
if [ $? -ne 0 ] ; then echo "ERROR: Exiting..." ; exit ; fi
echo "Testing puregauge results (if no output, everything OK!)"
for l in *.res ; do
  diff $l ../comp_puregauge_oldformat_40x9_b6.20.bin/$l -q
  if [ $? -ne 0 ] ; then
    echo "File $l differs (or is not there)!"
  fi
done
cd ..
