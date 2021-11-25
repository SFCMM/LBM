# V0.0.3
## Planned
- run basic 3D NS Case
- restart

# V0.0.2
## Planned
- run basic 2D NS Case
- load grid from disk
- Couette in 3D 
- Poiseuille flow in 3D
- add simple LPT solver
- capalilary simulation
- contact angle
- output analytical solution if wanted
- check in periodic boundaries if all cells have been linked
- get boundary const access to grid information
- pass all surfaces to the boundary manager
- store reference to the surfaces in each boundary
- allow comments in configuration files
- inflow pressure bc
- outflow pressure bc
- refactor forcing setup

# V0.0.1
## Planned:
- implement toRun()
- run basic 2D navier-stokes case (couette)
- load bcs from input file
- additional timers 
- switch analytical solution over to vector results!
- solve 2D poisson equation
- [bug] fix unused values being marked correctly!
- order boundary conditions by type
- exit on nan
- improve poiseuille setup (skip inlet and outlet cells for error calculation)

## Ongoing:
- implement output correctly


## Done:
### Features
- load grid generated from the generator
- setup diagonal neighbor connection for D2Q9 and similar methods
- periodic boundaries
- general bounce back "no-slip" wall boundary condition
- general bounce back "no-slip" wall boundary condition with tangential velocity
- load bnd conditions from configuration file
- added forcing

### Buildsystem

### Testing
- Added Couette flow testcase
- Added Poiseuille flow testcase

### Usability

### IO

### Performance

### Bugs

### Documentation

### Refactoring
