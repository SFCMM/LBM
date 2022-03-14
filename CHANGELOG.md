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
- [LPT] add simple collision model
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
- [LPT] calculate timestep by characteristic time
- [LBM][POISSON][BUG] the precision of D2Q5 is not correct for Poisson
- [LBM] change convergence criterium to some average not maximum

## moved

- load grid from disk
- Couette in 3D
- Poiseuille flow in 3D
- contact angle
- [LPT] add tutorial
- [LPT] compare solution at multiple time points
- determine analytical solutions for 2D poisson cases

# V0.0.2

## Planned

- output analytical solution if wanted
- check in periodic boundaries if all cells have been linked
- get boundary const access to grid information
- pass all surfaces to the boundary manager
- store reference to the surfaces in each boundary
- allow comments in configuration files
- refactor forcing setup
- refactor IO code
- make types definitions in include/common/constant
- move comparison with analytical solution to postprocessing
- fix keys in unused config value detection
- improve boundary condition with methods that use wet nodes
- make neem boundary condition available for NS
- add way to just iterate over leaf cells
- order boundary conditions by type

- add more openmp
- in cartesiangrid don't use the configuration.json but the accessor to mark the unused values correctly

- [BUG] D2Q5 still uses diagonal neighbors also this is not necessary
- simplify boundary by not making it necessary to copy the cellIds for each boundary separately
- generalize NEEM boundary condition
- fix poisson method to allow changes of relaxation
- reduce todos below 100

- warn when overwriting files
- give option to not overwrite files
- set output format for the solvers
- set output variables for the solvers
- check in output functions if fileformat ending is already included
- simplify writing output

- [LBM] add support for boundary ghost cells in periodic boundary
- add INFO_OUTPUT with physical time update

### this month:
- run basic 2D NS sphere case
- run basic 2D NS obstruction case

- [LPT] finish for release

- [BUG] All wetnode boundary condition case for poiseuille don't converge against a precise enough solution? (
  overprediction in the center??)


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

- Lattice Boltzmann Code
    - add D3Q19 and D3Q27
    - Boundary conditions
        - Wet node equilibrium method wall
        - Wet node non-equilibrium extrapolation method wall
        - Wet node non-equilibrium bounce-back method wall
        - BB Dirichlet bnd
        - AntiBB pressure bnd
        - Neumann NEEM bnd
    - analytical test:
        - Poisson:
            - 1D
            - 1D
        - Navier-Stokes:
            - Couette
            - Poiseuille

- Cartesian grid
    - allow generation of boundary ghost cells

- Poisson
    - 1D Reaction example added for 1D catalyst slab

### Buildsystem

### Testing

- check for divergence when debug is active
- write test report in CSV
- for analytical testcases allow alignment with the boundary

### Usability

- add pressure forcing term in periodic boundary condition

### IO

-set output dir and solution file name

### Performance

### Bugs
- segmentation fault after divergence
- write out full double precision
- [LBM][POISSON] fix D2Q9 for Poisson

### Documentation
- improve documentation of testcases

### Refactoring

- reduce the number of debug levels

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
