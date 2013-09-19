#!/bin/bash
mkdir vidtmp 
cd vidtmp
glretrace -b -s pic ../3dclusters.x.trace
find . -name "*.png" -type f -size -1000c -exec rm {} \;
x=0; for i in $(ls *.png); do counter=$(printf %05d $x); ln -s "$i" "$counter".PNG; x=$(($x+1)); done
avconv -f image2 -r 60 -i %05d.PNG -r 60 -codec:v libx264 -preset slow -crf 21 -y ../vid_out.mp4
cd ..
rm -r vidtmp
