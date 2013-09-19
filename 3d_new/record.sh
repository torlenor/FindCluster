#!/bin/bash
rm *.trace
apitrace trace ../bin/3dclusters.x -s 40 -n 11 3dcluster_40.list
#glretrace -s - 3dclusters.x.trace | ffmpeg -r 20 -f image2pipe -vcodec ppm -i pipe: -r 25 -qscale 31 -y output.ogg
#glretrace -s - 3dclusters.x.trace | avconv -r 80 -f image2pipe -c:v ppm -i pipe: -r 60 -qscale 1 -y output.mp4
#glretrace -s - 3dclusters.x.trace | avconv -r 80 -f image2pipe -vcodec ppm -i pipe: -r 60 -codec:v libx264 -y output.mp4
 
