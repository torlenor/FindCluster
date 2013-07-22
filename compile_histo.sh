#!/bin/bash

g++ -O3 -o bin/genpollhisto.x genpollhisto.cpp -lgsl -lgslcblas
#g++ -O3 -o bin/genpollhisto1d.x genpollhisto1d.cpp -lgsl -lgslcblas
g++ -O3 -o bin/genpollhisto1d_jack.x genpollhisto1d_jack.cpp -lgsl -lgslcblas
