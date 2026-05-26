# SwarmSearch

This repository contains the implementation and performance comparison of various **Particle Swarm Optimization (PSO)** algorithms developed for the AMSC project at Polimi.

### Implemented Algorithms

We implemented three different algorithms, namely

- **CHOPSO**
- **Simulated Annealing PSO**
- **ELPSO**

Each algorithm has been designed to support **arbitrary non-linear constraints**, ensuring found solutions adhere to specific problem requirements.

The three algorithms were packed in a fully-templated, header-only, dependency-less (with optional support for OpenMP) library. In order to run benchmarks, the [nlohmann/json](https://github.com/nlohmann/json) library must be available.

### Documentation

The full project documentation can be found in the [official documentation](https://etabeta1.github.io/particle-swarm-optimization/).

### Credits

This project was developed by [Mantoan Giorgio](https://github.com/GiorgioMantoan02), [Meshcheriakov Dmitrii](https://github.com/Dima765Me) and [Oggioni Andrea](https://github.com/etabeta1).

This work is an extension of a previous project originally developed by a team of five, including us three, and our former collaborators [Melzi Marco](https://github.com/marcomelzi) and [Salvi Paolo](https://github.com/paull194).
