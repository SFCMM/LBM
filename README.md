# LBM Solver

## Features
- Solve Navier-Stokes equation
- Solve Poisson-Boltzmann equation
- Lagrange particle tracking
   - drag
   - evaporation
   - injection


## License
- Code: [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
- Data/Documentation: [![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC_BY--NC_4.0-blue.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

## Requirements

* GNU 10+ or CLang 10+
* OpenMPI v4+
* HDF5 v1.12+
* CMake 3.20+

## Compile

1) create build dir
   ```mkdir build```
2) execute cmake
   ```cmake ..```
3) compile the binary
   ```make -j 6```
4) run bench
   ```./lbm --bench```
   
## Tests

Run all testcases:
* /test/run.sh

Clean testcase directories:
* /test/clean.sh

## Commandline options

### Debug (-d, --debug)

Set the debug level to run:

- 0 off (default)
- 1 on
- 2 all checks

For example to run with all debug checks:
```./gridgenerator --debug 2```

Determines how many extra checks and output is produced!

### Configuration File (-c, --config)

Set the path to the configuration that is to be used.

### Benchmark (-b, --bench)

Run a benchmark.

### Help (-h, --help)

### Version (-v, --version)

