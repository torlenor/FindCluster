# genres.sh v2.2.0 - 20130506 1123
#!/bin/bash

fraction="0"

Ns=40
Nt="2 3 4 5 6 7 8 9 10 11 12 14 16 18 20"

echo "This is genres.sh v2.0.0. Genering results files..." 

rm -f percc.res
echo "# Ns Nt percprob percproberr" >percc.res
for t in $Nt ; do cat ${Ns}x${t}/npercc_${Ns}x${t}_f${fraction}.res ; done >> percc.res

rm -f area.res
echo "# Nt area areaerr arealargestnonpercc arealargestnonperccerr largestclusterweight largestclusterweighterr       largestnonpercclusterweight largestnonpercclusterweighterr" > area.res
for t in $Nt ; do cat ${Ns}x${t}/area_${Ns}x${t}_f${fraction}.res | grep -v "^#" ; done >> area.res

rm -f cutrate.res_tmp cutrate.res
for t in $Nt ; do echo "$t `cat ${Ns}x${t}.log |  grep Cut | grep err`" >> cutrate.res_tmp ; done

echo "# Nt cut cuterr" > cutrate.res
cat cutrate.res_tmp | awk '{print $1 " " $4 " " $8}' >> cutrate.res

rm cutrate.res_tmp

# Cluster weight
rm -f clusterweight_tmp clusterweight.res clusterweight_err_tmp clusterweight_err_tmp2 clusterweight_tmp2
for t in $Nt ; do echo "$t `cat ${Ns}x${t}.log | grep Maximum | grep -v err`" >> clusterweight_tmp ; done
for t in $Nt ; do echo "$t `cat ${Ns}x${t}.log | grep Maximum | grep err`" >> clusterweight_err_tmp ; done

cat clusterweight_tmp | awk '{print $1 " " $13}' > clusterweight_tmp2
cat clusterweight_err_tmp | awk '{print $1 " " $13}' > clusterweight_err_tmp2

echo "# Nt clusterweight clusterweighterr" > clusterweight.res
join clusterweight_tmp2 clusterweight_err_tmp2 >>  clusterweight.res

rm clusterweight_tmp clusterweight_tmp2 clusterweight_err_tmp clusterweight_err_tmp2

# New cluster weight
echo "# Nt largestclusterweight largestclusterweighterr avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr largestnonpercclusterweight largestnonpercclusterweighterr" > clusterweight_new.res
for t in $Nt ; do cat ${Ns}x${t}/clustersize_${Ns}x${t}_f${fraction}.res | grep -v "^#" ; done >> clusterweight_new.res

# Polyakov loop results
rm -f poll.res
echo "# Nt poll pollerr" > poll.res
for t in $Nt ; do cat ${Ns}x${t}/poll_${Ns}x${t}_f${fraction}.res ; done >> poll.res

