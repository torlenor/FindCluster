#!/bin/bash

# if [ $# -lt 3 ] ; then
#         echo "Usage:"
#         echo "  ./start_all_cluster.sh startf deltaf endf"
#         exit
# fi

latdir=`basename \`pwd\``
echo latdir

for l in `ls -d T*` ; do
echo "Starting $l..."
	echo -e "\rRemoving old data..."
	rm -f $l/*.res
	rm -f $l/status.output
	rm -f $l/clrun
	rm -f $l/p*.o*
	rm -f $l/p*.e*

  	cd $l 
	echo -e "\tGenerating clrun file..."
	rm -f clrun
	touch clrun
	echo "#!/bin/csh" >> clrun
	echo "#$ -cwd" >> clrun
	echo "#$ -S /bin/csh" >> clrun
	echo "#$ -V" >> clrun
	# echo "#$ -q all.q,mpi" >> clrun
	echo "#$ -q all.q" >> clrun
	echo "#$ -l h_vmem=1.50G" >> clrun
	echo "#$ -N p${latdir}_${l}" >> clrun
	echo "setenv MV2_ENABLE_AFFINITY 0" >> clrun
	echo "cd `pwd`" >> clrun
	echo "./start_v_T.sh" >> clrun
	echo "exit 0" >> clrun
	# done with clrun file
	echo -e "\tSubmitting job..."
	qsub clrun
  cd ..
done
