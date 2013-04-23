#!/bin/bash
./genres.sh
rm fit.log
gnuplot plot_boxes_fit.gpl
gnuplot plot_boxes.gpl
gnuplot plot_cutrate.gpl
gnuplot plot_fractdim.gpl
gnuplot plot_percc.gpl

cat fit.log | grep -A4 -B2 "Final"
