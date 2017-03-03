# A simulation model for open-evolving systems

A simulation code for the model of open-evolving system.
[T. Shimada "A universal transition in the robustness of evolving open systems" Sci. Rep. 4: 4082 (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923212/).

# Usage

## Prerequisites

The main simulation program is written in standard C++.
To make a plot from the simulation result, Python3 and numpy and matplotlib libraries are needed.

## Compile & Run

To build the simulation code, just run

```
make
```

A binary program `eos.out` is made.

To run the program, set input parameters as command line arguments.
The following script executes the simulation and make plots.

```
./run.sh <model_TYPE> <#_of_interactions_per_species> <SimulationTime> <iseed>
```

## input parameters

- model_TYPE:
  - model_TYPE % 2 -> 0: fixed M, 1: flat degree distribution in (1,M)
  - (model_TYPE/2) % 2 -> 0: Gaussian link weight, 1: flat weight distribution
  - Therefore, 0: standard model, 1: flat degree, 2: flat link weight, 3: flat deg. & weight

For example,

```sh
./eos.out 0 10 10000 1234
```

## Output files

- timeseries.dat
- _output.json

# Installer

Run `install.sh` to compile the code and register it on OACIS.

```sh
git clone https://github.com/yohm/sim_eos_model.git && sim_eos_model/install.sh
```

# License

The MIT License (MIT)

Copyright (c) 2017 T. Shimada

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
