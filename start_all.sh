#!/bin/bash

for l in `ls -d f0.*` ; do
  echo "Starting $l..."
  cd $l
  screen -d -m ./start.sh
  cd ..
done
