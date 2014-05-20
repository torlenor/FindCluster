#!/bin/bash

# Percolation number vs. cut parameter
rm npercc.res
for l in `ls -d f0.* | sed s/^f//` ; do echo $l `cat f$l/npercc_80x20_f*.res | grep -v ^#` >> npercc.res  ; done

# Cluster radius vs. cut parameter
rm clusterradius.res
for l in `ls -d f0.* | sed s/^f//` ; do echo $l `cat f$l/clusterradius_80x20_f*.res | grep -v ^#` >> clusterradius.res  ; done

# Polyakov loop vs. cut parameter
rm poll.res
for l in `ls -d f0.* | sed s/^f//` ; do echo $l `cat f$l/poll_80x20_f*.res | grep -v ^#` >> poll.res  ; done
