#!/bin/bash
fraction=`basename \`pwd\` | sed s/^f//`
nconf=100

../findcluster.x -s 40 -t 10 -f $fraction -n $nconf ../pollinput.list | tee status.output
