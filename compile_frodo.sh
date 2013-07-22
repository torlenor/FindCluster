#!/bin/bash

icc -O3 -o bin/genpollhisto.x genpollhisto.cpp -I/usr/people/phyk/07schadh/local/include -L/usr/people/phyk/07schadh/local/lib -lgsl -L/software/Intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread 

icc -O3 -o bin/genpollev_from_config.x genpollev_from_config.cpp -I/usr/people/phyk/07schadh/local/include -L/usr/people/phyk/07schadh/local/lib -lgsl -L/software/Intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread
icc -O3 -o bin/genpollev_from_lpoll.x genpollev_from_lpoll.cpp -I/usr/people/phyk/07schadh/local/include -L/usr/people/phyk/07schadh/local/lib -lgsl -L/software/Intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread
icc -O3 -Wall -o bin/findcluster_debug.x findcluster.cpp -DDEBUG
icc -O3 -Wall -o bin/findcluster.x findcluster.cpp
