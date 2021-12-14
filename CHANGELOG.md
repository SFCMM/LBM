# V0.1
## Planned
- add python frontend
- add gpu acceleration

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

# V0.0.2

## Planned

- run basic 2D NS sphere case
- run basic 2D NS obstruction case
- load grid from disk
- Couette in 3D
- Poiseuille flow in 3D
- add simple LPT solver
- add simple injection model
- contact angle
- output analytical solution if wanted
- check in periodic boundaries if all cells have been linked
- get boundary const access to grid information
- pass all surfaces to the boundary manager
- store reference to the surfaces in each boundary
- allow comments in configuration files
- refactor forcing setup
- optionally cancel calculation when any macro value is nan
- refactor IO code
- make types definitions in include/common/constant
- set output separately for the solvers
- move comparison with analytical solution to postprocessing
- fix keys in unused config value detection
- improve boundary condition with methods that use wet nodes
- make neem boundary condition available for NS
- add way to just iterate over leaf cells
- order boundary conditions by type
- add documentation for the testcases
- simplify writing output
- check in output functions if fileformat ending is already included
- add more openmp
- in cartesiangrid don't use the configuration.json but the accessor to mark the unused values correctly
- move kdtree to common lib project.
- move line to common lib project.
- [BUG] D2Q5 still uses diagonal neighbors also this is not necessary
- add prefix names and output folder configuration options
- warn when overwritting files
- give option to not overwrite files
- [BUG] tangential velocity is not set correctly for D2Q5 for the couette case
- automatically run all the testcases and produce report
- simplify boundary by not making it necessary to copy the cellIds for each boundary separately
- generalize NEEM bounary condition
- fix poisson method to allow changes of relaxation
- allow setting write out precision for doubles for ASCII
- determine analytical solution for 2D poisson cases
- reduce to three debug levels

### priority:

- [BUG] segmentation fault after divergence
- inflow pressure bc
- outflow pressure bc

### moved:

- [BUG] fix Poisson D2Q9
- [BUG] fix D2Q5 poiseuille case
- use pressure boundary condition for poiseuille cases

## Done:

### Features

### Buildsystem

### Testing

### Usability

### IO

### Performance

### Bugs

### Documentation

### Refactoring

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

### Buildsystem

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

### Performance

### Bugs

### Documentation

### Refactoring
