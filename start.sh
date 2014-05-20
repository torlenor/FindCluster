#!/bin/bash
fraction=`basename \`pwd\` | sed s/^f//`
nconf=20

../findcluster.x -s 80 -t 20 -f $fraction -n $nconf ../spollinput.list | tee status.output
