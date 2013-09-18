#!/bin/bash
if [ "Darwin" = `uname` ]
then
	g++ 3dclusters.cpp -O3 -framework OpenGL -framework glut -o ./bin/3dclusters.x
else
	g++ 3dclusters.cpp -fopenmp -O3 -lglut -lGLU -o ./bin/3dclusters.x -lGL
fi
