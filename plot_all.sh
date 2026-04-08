#!/bin/bash

for file in output/snapshots/wf_*.dat
do
    echo "Plotting $file"
    gnuplot -c plot_density.gscript "$file"
done

echo "Plotting correlations and reaction probabilities"
gnuplot plot_corr.gscript
gnuplot plot_rxnprob.gscript
gnuplot -c plot_ref.gscript "output/psiref.dat"
gnuplot plot_pops.gscript
gnuplot plot_pot.gscript
echo "Done :)"