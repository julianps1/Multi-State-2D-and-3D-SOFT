# plot_theta_x4.gp
# Usage: gnuplot -e "dim='x';val=4.00787;tol=1e-6" plot_theta_x4.gp
# Default: slice in x at val

dim = "y"
val = 1.35294
#val = 0.023622
tol = 1e-6

set terminal pngcairo size 900,600 enhanced font 'DejaVuSans,12'
set output sprintf('output/pot_%s%.6g.png', dim, val)
set title sprintf('Diabatic Surfaces at %s=%.6g', dim, val)
set xlabel (dim eq "x") ? 'y' : 'x'
set ylabel 'Energy (E_h)'
set grid

# Data file columns: x y theta
# If slicing at constant x: filter by col1 and print y,theta ($2,$3)
# If slicing at constant y: filter by col2 and print x,theta ($1,$3)
awk_cond = (dim eq "x") ? "($1 - val) * ($1 - val) < tol*tol {print $2, $3}" : "($2 - val) * ($2 - val) < tol*tol {print $1, $3}"
awk_cond2 = (dim eq "x") ? "($1 - val) * ($1 - val) < tol*tol {print $2, $4}" : "($2 - val) * ($2 - val) < tol*tol {print $1, $4}"
plot sprintf("< awk -v val=%g -v tol=%g '%s' output/v_0.dat", val, tol, awk_cond) using 1:2 with l lw 3 lc "red" title "V00" , \
     sprintf("< awk -v val=%g -v tol=%g '%s' output/v_0.dat", val, tol, awk_cond2) using 1:2 with l lw 3 lc "blue" title "V11"

set output sprintf('output/adipot_%s%.6g.png', dim, val)
set title sprintf('Adiabatic Surfaces at %s=%.6g', dim, val)

# Data file columns: x y theta
# If slicing at constant x: filter by col1 and print y,theta ($2,$3)
# If slicing at constant y: filter by col2 and print x,theta ($1,$3)
awk_cond = (dim eq "x") ? "($1 - val) * ($1 - val) < tol*tol {print $2, $6}" : "($2 - val) * ($2 - val) < tol*tol {print $1, $6}"
awk_cond2 = (dim eq "x") ? "($1 - val) * ($1 - val) < tol*tol {print $2, $8}" : "($2 - val) * ($2 - val) < tol*tol {print $1, $8}"
plot sprintf("< awk -v val=%g -v tol=%g '%s' output/v_0.dat", val, tol, awk_cond) using 1:2 with l lw 3 lc "red" title "S0" , \
     sprintf("< awk -v val=%g -v tol=%g '%s' output/v_0.dat", val, tol, awk_cond2) using 1:2 with l lw 3 lc "blue" title "S1"


