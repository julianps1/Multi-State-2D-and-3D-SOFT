# qm2d_C

A small C++ code for propagating a two-state wavepacket on a 2D grid using FFT-based split-operator time evolution.

The project currently:

- reads simulation parameters from `INq`
- builds a fixed-size `256 x 256` spatial grid
- initializes a Gaussian wavepacket on one electronic state
- applies kinetic propagation with FFTW
- applies the potential step in a locally diagonalized adiabatic basis
- writes wavefunction densities, potential data, and correlation data to `output/`

## Current status

The codebase is actively in progress. 

So this README reflects the intended structure of the code as it exists now, plus the outputs and workflow already implemented.

## Requirements

- `g++` with C++17 support
- FFTW3 development libraries
- `make`
- optional: `gnuplot` for visualization

On a Linux system, the key dependency is the FFTW3 library used through `-lfftw3`.

## Build

From the project root:

```bash
make
```

This uses the provided `Makefile` and targets an executable named `sim`.

To remove object files and the executable:

```bash
make clean
```

## Run workflow

The intended workflow is:

1. Edit `INq` to set the grid, wavepacket, and propagation parameters.
2. Build with `make`.
3. Run `./sim`.
4. Inspect data files in `output/`.
5. Optionally generate plots with the supplied gnuplot scripts.

At startup, the code creates the `output/` directory if needed and deletes `output/correl.dat` so each run starts with a fresh correlation file.

## Input file

The default input file is `INq`. The parser reads whitespace-separated `name value` pairs.

The following parameters are currently recognized:

- `dt`: major time step
- `tmax`: number of outer time steps
- `tsub`: number of substeps per outer step
- `xmin`, `xmax`: x-grid bounds
- `ymin`, `ymax`: y-grid bounds
- `alfa`, `beta`: Gaussian width parameters
- `q0`, `y0`: initial wavepacket center
- `px0`, `py0`: initial momenta
- `nwpackets`: interval for writing wavefunction snapshots
- `h0`, `h1`, `h2`: mass-related parameters used to construct reduced masses
- `istate`: initial populated electronic state
- `pot_name`: potential selector parsed by the input reader, but not yet fully wired because `pick_pot(...)` is still missing

Current sample `INq`:

```text
dt 0.4
tmax 500
tsub 4
xmin -2.0
xmax 2.0
ymin -2.0
ymax 2.0
alfa 9.8
q0 1.0
px0 0
beta 9.8
y0 1.0
py0 0
nwpackets 100
h0 29
h1 1.0
h2 1.0
istate 0
```

## Output files

The code writes data files into `output/`.

### Wavefunction snapshots

Files:

- `output/wf_0.dat`
- `output/wf_1.dat`
- `output/wf_2.dat`
- ...

Each row contains:

- column 1: `x`
- column 2: `y`
- column 3: `|psi_0|^2`
- column 4: `|psi_1|^2`

Blank lines are inserted between x-slices for convenient `gnuplot pm3d` plotting.

### Correlation function

File:

- `output/correl.dat`

Each row contains:

- column 1: time `t`
- column 2: real part of the overlap on state 0
- column 3: imaginary part of the overlap on state 0
- column 4: real part of the overlap on state 1
- column 5: imaginary part of the overlap on state 1

This file is appended during propagation and reset at the beginning of each run.

### Potential data

File:

- `output/v_0.dat`

Current columns:

- column 1: `x`
- column 2: `y`
- column 3: diabatic `V00`
- column 4: diabatic `V11`
- column 5: diabatic coupling `V01`
- column 6: real part of lower adiabatic eigenvalue
- column 7: imaginary part of lower adiabatic eigenvalue
- column 8: real part of upper adiabatic eigenvalue
- column 9: imaginary part of upper adiabatic eigenvalue
- column 10: `cos(theta)`
- column 11: `sin(theta)`

### Mixing angle

File:

- `output/theta.dat`

Each row contains:

- column 1: `x`
- column 2: `y`
- column 3: local mixing angle `theta`

## Physics / numerical model

From the current implementation:

- the grid size is fixed at `ni = 256`, `nj = 256`
- the number of electronic states is fixed at `ns = 2`
- FFTW is used for 2D forward and backward transforms
- kinetic evolution is applied in momentum space
- the potential step is applied after diagonalizing the local `2 x 2` potential matrix
- propagation uses a symmetric split-operator pattern: half kinetic, full potential, half kinetic

The initial wavepacket is a Gaussian placed at `(q0, y0)` with widths controlled by `alfa` and `beta`, momentum shifts `px0` and `py0`, and initial occupation on state `istate`.

## Potential models

The source currently contains two potential-building functions:

- `pot_psi(GridData&)`
- `pot_Ferretti(GridData&)`

However, the selector function `pick_pot(...)` is not yet implemented in the current tree, so potential selection is not complete yet.

## Plotting

The repository includes simple gnuplot helpers:

- `plot_density.gscript`: plots state densities from `wf_*.dat`
- `plot_corr.gscript`: plots the correlation function
- `plot_pot.gscript`: plots potential-related columns from `v_*.dat`
- `plot_all.sh`: loops over `output/wf_*.dat` and generates density plots

Example intended usage:

```bash
gnuplot -c plot_corr.gscript
gnuplot -c plot_density.gscript output/wf_0.dat
gnuplot -c plot_pot.gscript output/v_0.dat
./plot_all.sh
```

## Source layout

- `src/main.cpp`: program entry point and time loop
- `src/globals.cpp`: input parsing from `INq`
- `src/init_wave.cpp`: reduced-mass setup
- `src/init_grid.cpp`: grid construction, initial wavefunction, potential diagonalization
- `src/FFT2D.cpp`: FFTW wrapper and transform execution
- `src/SO.cpp`: split-operator propagation
- `src/Ham.cpp`: initial energy expectation value helpers
- `src/Analysis.cpp`: output writers for densities, potentials, and correlations
- `src/pot_lib.cpp`: potential energy surface definitions
- `include/*.hpp`: declarations and global constants

## Notes for future cleanup

Natural next steps for the codebase are:

- define `grid::pot_name`
- implement `pick_pot(...)`
- document the meaning of each available `pot_name` value
- add a small example run that is known to compile and execute
- reduce reliance on global state as the code grows

