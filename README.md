# qm2d_C

Github repository: https://github.com/julianps1/Multi-State-2D-and-3D-SOFT

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
- Lapack development libraries
- `make`
- optional: `gnuplot` for visualization

On a Linux system, code requires the FFTW3 library used through `-lfftw3` and Lapack libraries `-llapacke -llapack -lblas`

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

0. Currently: set grid and number of states in globals.hpp, will need to `make clean` and recompile 
1. Edit `INq` to set the grid(later), wavepacket, and propagation parameters.
2. Build with `make`. (Only if constants in model potentials or grid were changed. In case of grid use make clean)
3. Run `./sim`.
4. Inspect data files in `output/`.
5. Plot :)

At startup, the code clears and creates the `output/` directories

## Input file

Outside of the input, in `globals.hpp` the number of gridpoints `ni,nj` and state `ns` are set globally (Make sure to do a full recompilation after modifying!!).

`ns` can be more than 2, however no potentials are currently defined with more than 2 states

The default input file is `INq`. The parser reads whitespace-separated `name value` pairs.

The following parameters are currently recognized:

- `dt`: major time step
- `tmax`: number of major time steps
- `tsub`: number of substeps per major step
- `xmin`, `xmax`: x-grid bounds
- `ymin`, `ymax`: y-grid bounds
- `alfa`, `beta`: Gaussian width parameters x,y
- `q0`, `y0`: initial wavepacket center x,y
- `px0`, `py0`: initial momenta x,y
- `nwpackets`: interval for writing wavefunction snapshots (snapshots every nwpackets steps)
- `h0`, `h1`, `h2`: mass-related parameters, for triatomic systems code uses jacobi coordinates, h0,h1,h2 correspond to atoms ABC for rxns A+BC,AB+C. If h2 is given as 0.0 code will use h0 and h1 as x and y masses  
- `istate`: initial populated electronic state
- `rCAPx,rCAPy`: Cutoffs for survival region and complex absorbing potential (CAO). If wanting a 1D cutoff, set the dimension you dont want a cutoff in to its max grid bound (i.e., if only interested in a cutoff in x, set rCAPy=ymax)
- `kCAP`: Strength of CAP, set to 0 for no CAP. IMPORTANT: SHOULD BE NEGATIVE!!
- `nmom`: Compute all moments up to `nmom` (can become very expensive). Set to 0 to skip
- `pot_name`: Choose your potential from `pot_lib.cpp`, 0-Landau Sombrero, 1-Ferretti, 2-Double Well, 3-NaFH, 4-H2+ model (Check that `ns` matches the chosen model, some only support 1 state!!)

Sample `INq`:

```text
dt 0.1 		
tmax 20000
tsub 1
xmin 1.0
xmax 8.0		
ymin 0.5
ymax 20.0
rCAPx 5.0
rCAPy 8.0
kCAP -2.0
alfa 4.0
q0 4.0
px0 0.0
beta 16.66
y0 1.76
py0 5.0
nwpackets 1000
h0 29.0
h1 19.0
h2 1.0
istate 1
nmom 0
pot_name 3
```
## Physics / numerical model

From the current implementation:

- FFTW is used for 2D forward and backward transforms
- kinetic evolution is applied in momentum space (and diabatic representation)
- the potential step is applied in the adiabatic representation (diagonalization and transformation automaticcaly accoutned for)
- Input potentials for the DIABATIC representation (no derivative coupling)
- propagation uses a symmetric split-operator pattern: half kinetic, full potential, half kinetic

The initial wavepacket is a Gaussian placed at `(q0, y0)` with widths controlled by `alfa` and `beta`, momentum shifts `px0` and `py0`, and initial occupation on state `istate`.

## Plotting

The repository includes simple gnuplot helpers in the `plotting` directory.

The `plot_all.sh` will execute all of them (regardless if there data files exist so some scripts may throw errors).

After plotting all wavefunction snapshots (with `plot_all.sh`), use `Make_gif.sh state` to create a gif of the state specified in the command.

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

- More dimensions
- reduce reliance on global state as the code grows
