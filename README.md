# Installation

1. Clone the repository

`git clone https://gitlab.rs.e-technik.tu-darmstadt.de/arctic/arctic.git`

2. Install the necessary Python modules.

<< List >>

3. Set the Python binary in **ARCTICsyn/sim.config**

`PYTHON_BINARY=<<Path to your Python binary>>`

# Run the First Synthesis with Standard Settings

The following example command calls ARCTIC via Gradle and synthesizes the three input Boolean function **a&(b|c)**.

`./ARCTICsyn/gradlew run --args="-f a&(b|c)"`

This command synthesizes the circuit structures and performs technology mapping on them with the preset genetic gate library. When the process is finished, the best found circuit structure and assignment of logic gates as well as the corresponding score is printed out. In the directory **ARCTICsyn/benchmarks** a new directory corresponding to the run has been created containing all found circuit structures as DOT and JSON files as well as structure and assignment of gates for the resulting circuit.

# Setting for Synthesis, Simulation and Technology Mapping

## Synthesis

The synthesis settings are located in the configuration file **ARCTICsyn/syn.config**.

### OUTPUT_DIR

### SYN_LIMIT_THREAD_NUM

### SYNTHESIS_DEPTH

### SYNTHESIS_WEIGHT

### SYNTHESIS_WEIGHT_RELAXATION

### SYNTHESIS_LIMIT_STRUCTURES_NUM

## Simulation

The simulation settings are located in the configuration file **ARCTICsyn/sim.config**.

### PYTHON_BINARY
### SIM_LIMIT_THREADS_NUM
### SIM_INIT_ARGS
### SIM_ARGS

## Technology Mapping

The technology mapping settings are located in the configuration file **ARCTICsyn/map.config**.

### LIBRARY
### COMPAT_LIBRARY
### SEARCH_ALGORITHM
### OPTIMIZATION_TYPE
### STATISTICS
### BAB-SEARCH_STRATEGY
### BAB-TYPE
### BAB-VISUALIZE
### BAB-STATISTICS
### BAB-FAST
