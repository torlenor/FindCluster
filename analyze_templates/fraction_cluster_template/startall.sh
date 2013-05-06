#!/bin/bash
echo "Starting everything!"
read -p "Are you sure? " -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
 for f in `ls -d f0.*` ; do cd $f ; rm -r 40x* analyze* c40x*.sh ; ./genp.sh ; for l in `ls c40x*` ; do qsub $l ; done ; cd .. ; done
fi
