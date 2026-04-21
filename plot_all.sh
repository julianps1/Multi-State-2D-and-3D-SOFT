#!/bin/bash

for file in output/snapshots/wf_*.dat
do
    echo "Plotting $file"
    gnuplot -c plotting/plot_density.gscript "$file"
done

echo "Plotting correlations and reaction probabilities"
gnuplot plotting/plot_corr.gscript
gnuplot plotting/plot_rxnprob.gscript
gnuplot -c plotting/plot_ref.gscript "output/psiref.dat"
gnuplot plotting/plot_pops.gscript
gnuplot plotting/plot_pot.gscript
gnuplot plotting/plot_moments.gscript
echo "Done :)"