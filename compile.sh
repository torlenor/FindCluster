#!/bin/bash

# g++ -O3 -o bin/genpollev_from_config.x genpollev_from_config.cpp -lgsl -lgslcblas -llapack
# g++ -O3 -o bin/genpollev_from_lpoll.x genpollev_from_lpoll.cpp -lgsl -lgslcblas -llapack
# g++ -O3 -Wall -o bin/findcluster_debug.x findcluster.cpp -DDEBUG
g++ -g -O3 -Wall -o bin/findcluster.x findcluster.cpp
