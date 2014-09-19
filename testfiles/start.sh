#/bin/bash
/usr/bin/time ./findcluster.x -s 32 -t 8 -n 5 qcd.list -b -a -w -f 0.00 > qcd_status_f0.00.out
/usr/bin/time ./findcluster.x -s 32 -t 8 -n 5 qcd.list -b -a -w -f 0.45 > qcd_status_f0.45.out
mkdir tmp_qcd
mv *.res tmp_qcd
mv qcd_status_f0.*.out tmp_qcd
/usr/bin/time ./findcluster.x -s 40 -t 9 -n 5 puregauge.list -b -a -w -f 0.00 -o > puregauge_status_f0.00.out
/usr/bin/time ./findcluster.x -s 40 -t 9 -n 5 puregauge.list -b -a -w -f 0.45 -o > puregauge_status_f0.45.out
mkdir tmp_puregauge
mv *.res tmp_puregauge
mv puregauge_status_f0.*.out tmp_puregauge
