# V0.1

## Planned

- add python frontend
- add gpu acceleration
- reduce todo to 0
- [LPT] eliptic particles
- [LPT] fluid collision model
- [LPT] breakup model
- [LPT] dynamic drag model
- [LPT] plastic collision

# V0.0.4

## Planned

- [LPT] calculate timestep by characteristic time
- [LPT] add simple collision model

# V0.0.3

## Planned

- run basic 3D NS Case
- restart
- store meta data in output files (version, node etc.)
- output performance data
- capillary simulation
- add support for multiple species
- add support for multiple phases
- add injection models
- add breakup models
- add collision code for LPT
- add multilevel method
- add github website
- add images to readme
- allow more than one solver to run at the same time
- write out lines as svg
- set write out precision
- write out memory statistics
- reduce todo below 50

- ray tracing of particles
- refactor comparison with analytical solution to postprocessing
- initial material point method
- mpi communication methods
- [LBM] calculate total energy/internal energy for output purposes (also output pressure and temperature)
- compare timings between runs of the testing tool
- port testing tool over to python
- [LIB] move more stuff to common lib project
- [LIB] refactor kdtree and add more complete test
- [LBM][POISSON][BUG] the precision of D2Q5 is not correct for Poisson
- [LBM] change convergence criterium to some average not maximum
- [LBM] fix NEBB bnd cnd
- [LBM] add zou-he bnd cnd
- [LBM][BUG] Periodic bnd doesn't work always (see step_ns_bug1.json)
- Allow in IO functions to set output variable types
- [test script] in test script check if bin is current
- [test script] run multiple bins 
- [test script] backup binary for each run
- [test script] run test on remote host
- update to gnu 11
- update to clang 11
- add intel compiler option

## moved

- load grid from disk
- Couette in 3D
- Poiseuille flow in 3D
- contact angle
- [LPT] add tutorial
- [LPT] compare solution at multiple time points
- determine analytical solutions for 2D poisson cases
- [BUG] All wetnode boundary condition case for poiseuille don't converge against a precise enough solution?
  (overprediction in the center??) (related to the uncertainty about calculating the relaxation??)
- move comparison with analytical solution to postprocessing
- order boundary conditions by type (is this needed?)
- [BUG] D2Q5 still uses diagonal neighbors also this is not necessary(is this needed?)
- [LBM] add support for boundary ghost cells in periodic boundary (needed??)
- [IO] allow switching of output filetypes by handling a configuration object
- [IO] simplify writing output
- [IO] warn when overwriting files
- [IO] give option to not overwrite files
- [IO] set output format for the solvers
- [IO] set output variables for the solvers
- [IO] check in output functions if fileformat ending is already included
- add more openmp
- allow comments in configuration files
- add ns compatible neumann bndcnd
- test all bnds for both poisson and NS
    - working:

# V0.0.2

## Planned

- output analytical solution if wanted
- get boundary const access to grid information
- pass all surfaces to the boundary manager
- store reference to the surfaces in each boundary
- refactor forcing setup
- fix keys in unused config value detection
- add way to just iterate over leaf cells
- in cartesiangrid don't use configuration.json but the accessor to mark the unused values correctly
- reduce todos below 100

### Working on:
- make types definitions in include/common/constant


## Done:

### Features

- Lagrange Particle Solver
    - Forces:
        - Gravity
        - Buoyancy
        - Stokes Drag
        - Non-linear Drag
    - Integration methods:
        - Forward Euler
        - Implicit Euler
    - analytical tests:
        - free fall
        - stokes drag
        - terminal velocity
    - simple injection model

- Lattice Boltzmann Solver
    - add D3Q19 and D3Q27
    - Boundary conditions
        - Wet node equilibrium method wall
        - Wet node non-equilibrium extrapolation method wall
        - Wet node non-equilibrium bounce-back method wall
        - BB Dirichlet bnd
        - AntiBB pressure bnd
        - Neumann NEEM bnd
    - Analytical test:
        - Poisson:
            - 1D Debye-Hueckel
            - 1D Diffusion Slab
            - 2D Helmholtz Equation
        - Navier-Stokes:
            - Couette
            - Poiseuille
  - Examples:
      - Navier-Stokes:
          - single step
          - double step
          - flow around cylinder
      - Poisson
          - single step
          - double step
          - 2D diffusion
          - 2D Debye-Hueckel

- Math expression support for boundary conditions

### Buildsystem

### Testing
- check for divergence when debug is active
- write test report in CSV
- for analytical testcases allow alignment with the boundary

### Usability
- add pressure forcing term in periodic boundary condition
- allow to set a bnd per geometry obj
- allow setting "conv_interval" for the interval for the convergence check
- write output every "60" seconds by default as a keep alive msg
    - time interval is settable by "keep_alive_time"

### IO
- set output dir and solution file name
- allow multiple cell output filters
- correctly handle switch between int and double in writePoints()

### Performance

### Bugs
- segmentation fault after divergence
- write out full double precision
- [LBM][POISSON] fix D2Q9 for Poisson
- error in Base64 encoding

### Documentation
- improve documentation of testcases

### Refactoring

- reduce the number of debug levels
- refactor output cell filter code
- remove shift from base64 which didn't work

# V0.0.1

### Features

- load grid generated from the generator
- setup diagonal neighbor connection for D2Q9 and similar methods
- periodic boundaries
- general bounce back "no-slip" wall boundary condition
- general bounce back "no-slip" wall boundary condition with tangential velocity
- Poisson Non-equilibrium extrapolation method based dirichlet boundary condition
- load bnd conditions from configuration file
- added forcing
- Postprocessing
    - Lines
- basic openMP
- supports D1Q3, D2Q5 and D2Q9
- Support the solution of Navier-Stokes and Poisson-Boltzmann equation

### Testing

- Added Couette flow testcase
- Added Poiseuille flow testcase
- option to exclude surfaces from error calculation

### Usability

- exit on nan during convergence check
- Add postprocessing structures
- define boundary condition in the configuration file

### IO

- differentiate between float and integers correctly in the VTK output
- write lines output to csv
