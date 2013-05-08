#!/bin/bash

dir=/usr/people/phyk/07schadh/cltmp/cluster_20130506/analyze/Ns40/b6.20

rsync --exclude=*.gpl --exclude=analyze* --exclude=*.x -avz frodo:${dir}/* .
