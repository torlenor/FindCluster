#!/bin/bash

icc -O3 -o bin/genpollhisto.x genpollhisto.cpp -I/usr/people/phyk/07schadh/local/include -L/usr/people/phyk/07schadh/local/lib -lgsl -L/software/Intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread 
icc -O3 -o bin/gen3dpolllattice.x gen3dpolllattice.cpp -I/usr/people/phyk/07schadh/local/include -L/usr/people/phyk/07schadh/local/lib -lgsl -L/software/Intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread
