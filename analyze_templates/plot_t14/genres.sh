#!/bin/bash

rm -f percc.res
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do cat 40x${t}/npercc_40x${t}.res ; done > percc.res

rm -f cutrate.res_tmp cutrate.res
for t in 2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 ; do echo "$t `cat 40x${t}.log |  grep Cut | grep err`" >> cutrate.res_tmp ; done

echo "# Nt cut cuterr" > cutrate.res
cat cutrate.res_tmp | awk '{print $1 " " $4 " " $8}' >> cutrate.res

rm cutrate.res_tmp
