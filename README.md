# SwarmSearch
<!-- Optimization by swarm search -->

## Description of the project
The aim of the project is to implement the Particle Swarm Optimization (PSO) method

Particle swarm optimization is a computational method 
that mimics the movement of a flock of birds or a swarm of flies.
* The method uses a population of candidate solutions, which are dubbed particles.
* These particles move around in the search-space according to simple laws based on their position.
* Each particle's movement is influenced by its local best known position (personal best solution) and is also guided toward the best known positions in the search-space (global best solution), which are updated as better positions are found by other particles.

In practice we aim to find the global best solution in this way:

$global_{best} = {argmin}_{x \in D} f(x)$, 
where D is a subspace of $R^n$

## Goals of the project
* Implement a serial version of one of the variants of the code.
* Test the code on a set of hard minimization problems.
* Parallelize the code using OpenMP or MPI.

## The Algorithm
PSO optimizes a problem by iteratively improving a candidate solution with regard to a given measure of quality.

### Mathematical Formulation
The swarm consists of $N$ particles. At each iteration $t$, every particle $i$ updates its:
1.  **Velocity ($v$):** Based on inertia, personal position, and swarm position.
2.  **Position ($x$):** Based on the new velocity.

The update equations are defined as:

**1. Velocity Update:**
$v_{i}^{t+1} = \omega v_{i}^{t} + c_{1}r_{1}(p_{best_{i}}^{t}-x_{i}^{t}) + c_{2}r_{2}(g_{best}^{t}-x_{i}^{t})$

**2. Position Update:**
$x_{i}^{t+1} = x_{i}^{t} + v_{i}^{t+1}$

**Where:**
* $v_{i}^{t}$: Current velocity of particle $i.
* $x_{i}^{t}$: Current position of particle $i$.
* $\omega$: Inertia weight.
* $p_{best}$: The best position found by the particle so far (Personal Best Solution).
* $g_{best}$: The best position found by the entire swarm so far (Global Best Solution).
* $c_1, c_2$: Acceleration coefficients (constants).
* $r_1, r_2$: Random vectors uniformly distributed in $[0,1]$.

### Pseudocode
The logic implemented in this project follows this structure:

```
Initialize Swarm:
  For each particle i in N:
    Initialize position x_i and velocity v_i randomly
    Calculate fitness f(x_i)
    Set p_best_i = x_i
    if f(x_i) < g_best
        g_best = x_i

Loop until Max Iterations (T) or Convergence:
  For each particle i in N:
    // Update Velocity
    v_i = w*v_i + c1*r1*(p_best_i - x_i) + c2*r2*(g_best - x_i)
    
    // Update Position
    x_i = x_i + v_i
    
    // Evaluation
    If f(x_i) < f(p_best_i):
      p_best_i = x_i
      If f(p_best_i) < f(g_best):
        g_best = p_best_i
END
