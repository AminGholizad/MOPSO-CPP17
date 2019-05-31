# Multi-Objective Particle swarm optimization (MOPSO) CPP17
Multi-Objective particle swarm optimization algorithm ([MOPSO](https://en.wikipedia.org/wiki/Particle_swarm_optimization#Multi-objective_optimization)) for a minimization problem. In this project, nonlinar constraints are implemented as infeasable solutions.
This project is implemented in C++17.
Constrtaint and objectives are in one function.

# Features

1. Template class for particles with number of variables, objectivess and swarm size as the template parameters is used to initialize the arrays of the needed size.
2. Template class for repository with number of variables,objectives and swarm size as the template parameters is used to initialize the arrays of the needed size.
3. Template function for pso with number of variables and swarm size is as the template parameters is used to initialize arrays accordingly
4. The cost function should accept one array as input and output a pair of doubles representing cost and constaint of the problem respectively.
5. The Mutiation is used to avoid local minima.

# todo
[] More explanation of the functionality
