# LBM Solver

## Features
- Solve Navier-Stokes equation
- Solve Poisson-Boltzmann equation
- Lagrange particle tracking
   - drag
   - evaporation
   - injection


## License
- Code is BSD v3 licensed
- Data and Documentation is CC-BY-NC ![license](https://creativecommons.org/licenses/by-nc/4.0/ "CC-BY-NC)")

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

## Commandline options

### Debug (-d, --debug)

Set the debug level to run:

- 0 off (default)
- 1 minimal
- 2 normal
- 3 more
- 4 maximum

For example to run with the debug level 3:
```./gridgenerator --debug 3```

Determines how many extra checks and output is produced!

### Configuration File (-c, --config)

Set the path to the configuration that is to be used.

### Benchmark (-b, --bench)

Run a benchmark.

### Help (-h, --help)

### Version (-v, --version)

