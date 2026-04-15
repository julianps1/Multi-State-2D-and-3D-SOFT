# plot_theta_x4.gp
# Usage: gnuplot -e "dim='x';val=4.00787;tol=1e-6" plot_theta_x4.gp
# Default: slice in x at val

dim = "y"
val = 0.023622
tol = 1e-6

set terminal pngcairo size 900,600 enhanced font 'DejaVuSans,12'
set output sprintf('output/theta_%s%.6g.png', dim, val)
set title sprintf('Theta slice at %s=%.6g', dim, val)
set xlabel sprintf('%s', dim)
set ylabel 'theta (rad)'
set grid

# Data file columns: x y theta
# If slicing at constant x: filter by col1 and print y,theta ($2,$3)
# If slicing at constant y: filter by col2 and print x,theta ($1,$3)
awk_cond = (dim eq "x") ? "($1 - val) * ($1 - val) < tol*tol {print $2, $3}" : "($2 - val) * ($2 - val) < tol*tol {print $1, $3}"
plot sprintf("< awk -v val=%g -v tol=%g '%s' output/theta.dat", val, tol, awk_cond) using 1:2 with linespoints lw 1 pt 7 title sprintf('%s=%.6g', dim, val)
