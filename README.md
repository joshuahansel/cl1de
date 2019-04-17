# Introduction

`cl1de` stands for "**C**onservation **L**aw **1D** **E**xplicit".
As its name suggests, it is used to solve conservation laws in 1D,
using explicit time integration. Currently the following
conservation laws are supported:
* Variable-area compressible, inviscid Euler equations

The following conservation laws are also planned:
* 2-phase, variable-area, compressible, inviscid Euler equations
* Shallow water equations

# Numerical Methods

Planned spatial discretization schemes include various flavors
of the following families:
* Finite volume
* Discontinuous Galerkin finite element

Temporal discretizations include the following:
* Explicit TVD Runge-Kutta (AKA SSPRK) methods of orders 1, 2, and 3

# Language

It is written in C++ with minimal external dependencies.

# Compilation

Compilation is performed using `make`. An executable is created
per problem:
```
make <problem_name>
```
where `<problem_name>` is the name of the problem. Each problem
resides in a first-level sub-directory of the `problems` directory,
and the associated main C++ file is named as `<problem_name.cpp>`.

# Execution

`make <problem_name>` creates the executable `problems/<problem_name>/run`.
The executable is run with the name of an input file as an argument:
```
./run <input_file>
```
The output is the file `output.csv`, which contains the final solution
vectors for the associated conservation law.
