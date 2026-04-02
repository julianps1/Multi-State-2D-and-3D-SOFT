# qm2d_C

A small C++ code for propagating a two-state wavepacket on a 2D grid using FFT-based split-operator time evolution.

The project currently:

- reads simulation parameters from `INq`
- initializes a Gaussian wavepacket on one electronic state
- applies kinetic propagation with FFTW
- applies the potential step in a locally diagonalized adiabatic basis
- writes wavefunction densities, potential data, and correlation data to `output/`
- works for a double well coupled to a harmonic oscillator 
## Current status

The codebase is actively in progress. 

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

0. Currently: set grid and number of states in globals.hpp
1. Edit `INq` to set the grid(later), wavepacket, and propagation parameters.
2. Build with `make`. (Only if constants in model potentials or grid were changed. In case of grid use make clean)
3. Run `./sim`.
4. Inspect data files in `output/`.
5. Plot :)

At startup, the code clears and creates the `output/` directories

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

Sample `INq`:

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
pot_name 2
```

## Output files

The code writes data files into `output/`, this is cleared before each run so SAVE YOUR DATA EXTERNALLY

### Wavefunction snapshots

Files:

- `output/snapshots/wf_0.dat`
- `output/snapshots/wf_1.dat`
- `output/snapshots/wf_2.dat`
- ...

Each row contains:

- column 1: `x`
- column 2: `y`
- column 3: `|psi_0|^2`
- column 4: `|psi_1|^2`

Blank lines are inserted between x-slices for convenient `gnuplot pm3d` plotting.

The total number is controlled by nwpackets, which sets the number of time-steps between each snapshot

### Correlation function (Auto and Cross)

File:

- `output/correl.dat` or `output\crosscorrel.dat`

Each row contains:

- column 1: time `t`
- column 2: real part of the overlap on state 0
- column 3: imaginary part of the overlap on state 0
- column 4: absolute value

The difference is the reference function. Autocorrelation is always psi(0) but cross correlation is determined by psiref in `init_psiref`
This is MODEL SPECIFIC

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

This is for debugging

File:

- `output/theta.dat`

Each row contains:

- column 1: `x`
- column 2: `y`
- column 3: local mixing angle `theta`

## Physics / numerical model

From the current implementation:

- the grid size is fixed and set in `globals.hpp` along with the number of states
- the number of electronic states is limited to 2 maximum
- FFTW is used for 2D forward and backward transforms
- kinetic evolution is applied in momentum space (and diabatic representation)
- the potential step is applied in the adiabatic representation (diagonalization and transformation automaticcaly accoutned for)
- Input potentials for the DIABATIC representation (no derivative coupling)
- propagation uses a symmetric split-operator pattern: half kinetic, full potential, half kinetic

The initial wavepacket is a Gaussian placed at `(q0, y0)` with widths controlled by `alfa` and `beta`, momentum shifts `px0` and `py0`, and initial occupation on state `istate`.

## Potential models

- 0 `pot_psi(GridData&)`
- 1 `pot_Ferretti(GridData&)`
- 2 `pot_doublewell(GridData&)` 

The values here correspond to the `pot_name` in the input. NOTE `pot_doublewell` is a 1-state model. Currently switching between 1 and two states is not automatic. 

## Plotting

The repository includes simple gnuplot helpers:

- `plot_density.gscript`: plots state densities from `wf_*.dat`
- `plot_corr.gscript`: plots the correlation function
- `plot_pot.gscript`: plots potential-related columns from `v_*.dat`
- `plot_all.sh`: loops over `output/wf_*.dat` and generates density plots and .gif

Example intended usage:

```bash
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

- Automatic switching between 1 and 2 states
- More potentials
- More dimensions
- add a small example run that is known to compile and execute
- reduce reliance on global state as the code grows
