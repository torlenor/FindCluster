#!/bin/bash
Ns=32
Nt=8
fraction=0.71367188
nconf=200

../findcluster.x -s $Ns -t $Nt -f $fraction -n $nconf ./pollinput.list -b | tee status.output
