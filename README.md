# A simulation model for open-evolving systems

A simulation code for the model of open-evolving system.
[T. Shimada "A universal transition in the robustness of evolving open systems" Sci. Rep. 4: 4082 (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3923212/).

# Compile & Run

Run `make`. `eos.out` is made.

To run the program, set input parameters as command line arguments.

```
./eos.out <model_TYPE> <#_of_interactions_per_species> <SimulationTime> <iseed>
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

Run `install.sh` will compile the codes and then registers the simulator on OACIS.

```sh
git clone https://github.com/yohm/sim_eos_model.git && sim_eos_model/install.sh
```
