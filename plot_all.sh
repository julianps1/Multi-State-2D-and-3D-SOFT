#!/bin/bash

for file in output/wf_*.dat
do
    echo "Plotting $file"
    gnuplot -c plot_density.gscript "$file"
done
