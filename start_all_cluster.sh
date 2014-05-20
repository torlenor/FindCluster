#!/bin/bash

if [ $# -lt 3 ] ; then
        echo "Usage:"
        echo "  ./start_all_cluster.sh startf deltaf endf"
        exit
fi

start=$1
delta=$2
end=$3

for f in `seq -f %.2f $start $delta $end` ; do
echo "Starting $f..."
  if [ ! -d f$f ] ; then
	mkdir f$f
  else
	rm -f f$f/*.res
	rm -f f$f/status.output
	rm -f f$f/clrun
  fi
  cd f$f
	ln -s -f  ../start.sh .
	ln -s -f ../spollinput.list .
	echo -e "\tGenerating clrun file..."
	rm -f clrun
	touch clrun
	echo "#!/bin/csh" >> clrun
	echo "#$ -cwd" >> clrun
	echo "#$ -S /bin/csh" >> clrun
	echo "#$ -V" >> clrun
	echo "#$ -q all.q,mpi" >> clrun
	echo "#$ -l h_vmem=3.75G" >> clrun
	echo "#$ -N spc${f}_8020" >> clrun
	echo "cd `pwd`" >> clrun
	echo "./start.sh" >> clrun
	echo "exit 0" >> clrun
	# done with clrun file
	echo -e "\tSubmitting job..."
	qsub clrun
  cd ..
done
