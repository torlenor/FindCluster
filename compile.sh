#!/bin/bash

g++ -O3 -o bin/genpollhisto.x genpollhisto.cpp -lgsl -lblas
g++ -O3 -o bin/gen3dpolllattice.x gen3dpolllattice.cpp -lgsl -lblas -llapack
