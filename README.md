# Installation

1. Clone the repository

`git clone https://gitlab.rs.e-technik.tu-darmstadt.de/arctic/arctic.git`

2. Install the necessary Python modules.

- numpy
- scipy
- autograd
- matplotlib


`pip install numpy scipy autograd matplotlib `

3. Set the Python binary in **ARCTICsyn/sim.config**

`PYTHON_BINARY=<<Path to your Python binary>>`

# Run the First Synthesis with Standard Settings

The following example command calls ARCTIC via Gradle and synthesizes the three input Boolean function **a&(b|c)**.

`./ARCTICsyn/gradlew run --args="-f a&(b|c)"`

This command synthesizes the circuit structures and performs technology mapping on them with the preset genetic gate library. When the process is finished, the best found circuit structure and assignment of logic gates as well as the corresponding score is printed out. In the directory **ARCTICsyn/benchmarks** a new directory corresponding to the run has been created containing all found circuit structures as DOT and JSON files as well as structure and assignment of gates for the resulting circuit.

# Setting for Synthesis, Simulation and Technology Mapping

## Synthesis

The synthesis settings are located in the configuration file **ARCTICsyn/syn.config**.

**OUTPUT_DIR** Output directory for synthesis results relative to ./ARCTICsyn/.

**SYN_LIMIT_THREAD_NUM** Limit for the number of threads used by the synthesis. When set to 0, the number of processor cores - 1 is used.

**SYNTHESIS_DEPTH** Maximum depth of the synthesized circuits.

**SYNTHESIS_WEIGHT** Maximum weight (number of gates) of the synthesized circuits.

**SYNTHESIS_WEIGHT_RELAXATION** Limit for the circuit weight w.r.t. the minimum circuit size. E.g. a weight relaxation of 2 allows circuits to be at most two gates bigger than the smallest circuit found.

**SYNTHESIS_LIMIT_STRUCTURES_NUM** Absolute limit for the generated number of circuits.

## Simulation

The simulation settings are located in the configuration file **ARCTICsyn/sim.config**.

**PYTHON_BINARY** Path to the Python binary to be used.

**SIM_LIMIT_THREADS_NUM** Limit for the number of threads used by the synthesis. When set to 0, the number of processor cores - 1 is used.

**SIM_INIT_ARGS** Initialization arguments for the simulator.

**SIM_ARGS** Simulation specific arguments for the simulator.

## Technology Mapping

The technology mapping settings are located in the configuration file **ARCTICsyn/map.config**.

**LIBRARY** Library of genetic gates to be used.

**COMPAT_LIBRARY** Compatibility library to be used.

**SEARCH_ALGORITHM** Search algorithm to be used. Currently supported: **EXHAUSTIVE, ANNEALING, BRANCH_AND_BOUND**

**STATISTICS** If true, statistics are generated and output in JSON format.

**BAB-SEARCH_STRATEGY** *(Branch-and-Bound specific)* Search strategy of B&B. Currently supported: **DEPTH_FIRST_SEARCH, BREADTH_FIRST_SEARCH, CYCLIC_BEST_FIRST_SEARCH, BEST_FIRST_SEARCH**

**BAB-TYPE** *(Branch-and-Bound specific)* Type of B&B.  Currently supported: **EAGER, LAZY**

**BAB-VISUALIZE** *(Branch-and-Bound specific)* If true, the search tree is visualized as DOT file.

**BAB-STATISTICS** *(Branch-and-Bound specific)* If true, B&B specific statistics are generated.

**BAB-FAST** *(Branch-and-Bound specific)* If true, heuristic mode is activated.

# Genetic Gate Libraries

Possible choices for the argument **LIBRARY** in **ARCTICsyn/map.config** are located in the folder **ARCTICsim/thermo_libs**. The gate libraries contain thermodynamic parameters for promoters, transcription factors (TFs) and the host context. With these parameters, circuit simulations can be performed. The core of the set of libraries is the 'ideal' or base library. It contains numerically fitted thermodynamic parameters for the 40 different gates available in Cello's library (`https://github.com/CIDARLAB/cello`).

## Library Generation

The different gate libraries are generated from the base library using different hyperparameters governing the distribution of crosstalk. The naming convention is the following.

`thermo_lib_id_xtalk_all_dirichlet_<DIRICHLET>tot<TOTAL>.json`

where
- `<DIRICHLET>` is a label of the concentration parameter of the Dirichlet distribution that governs the crosstalk distribution across non-cognate transcription factors in a gate. It is a value in `{0, 1, ..., 6}`. See below for the actual concentration values.
- `<TOTAL>` is a value specifying how strong crosstalking TF's modify the gate's output in relation to the cognate TF. It is a value in `[0, 1]`, where a `1` specifies that if all crosstalking TF's are present in a unit concentration, the gate's output would be modified as if there was instead the cognate TF present in a unit concentration.

The distribution of crosstalk is governed by a Dirichlet distribution. Thus, the distribution has support on all unit-length vectors with as many entries as there are non-cognate TFs. Usually, every entry has an individual concentration parameter assigned. In our case, these are all equal to a global concentration parameter. This specific instance of Dirichlet distribution is also often called a 'symmetric Dirichlet distribution'. The labels `{0, 1, ..., 6}` of this global concentration are mapped to real concentration values by the following table

label | real concentration (rounded)
--------|-----------------------------
`0`     |   0.004717
`1`     |   0.114532
`2`     |   0.378811
`3`     |   1.
`4`     |   2.63984
`5`     |   8.731212
`6`     | 212.005864

The lowest label `0` assigns most crosstalk to only one non-cognate TF while the highest label `6` gives a mostly uniform distribution of crosstalk across non-cognate TFs.
