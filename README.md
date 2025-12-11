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

## The Standard Algorithm
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
```

## Chaotic Particle Swarm Variant

### Overview
This project explores a variant of the standard PSO called **Chaotic Particle Swarm Optimization**. 

In this algorithm, we introduce, along with standard particles, chaotic particles which move in the dominion of the problem using chaotic maps, which allow to better explore all the dominion of the problem. Chaos is characterized by two distinctive properties: **non-repeatability** and **periodicity**.


### Main Differences: Standard vs. Chaotic

The key differences and improvements of the Chaotic variant over the Standard PSO are:

| Feature | Standard PSO | Chaotic PSO |
| :--- | :--- | :--- |
| **Particles** | Just normal particles. | Chaotic and normal particles. |
| **Search Speed** | Slower general search speed compared to chaotic methods. | Has the ability to perform general searches at **higher speeds**. |
| **Convergence** | Often haunted by **fruitless early convergence**. | The introduction of chaotic particles **avoid early convergence** through adaptability and ideal global search capabilities. |

### Modifications to the standard Algorithm
Compared to the standard algorithm, in the chaotic variant we have to distinguish between normal and chaotic particles, this is done in practice using two different sub-classes but we could modify the pseudocode in this way:

```
Loop until Max Iterations (T) or Convergence:
  For each particle i in N:
   if (i is normal)
   {
    // Update Velocity
    v_i = w*v_i + c1*r1*(p_best_i - x_i) + c2*r2*(g_best - x_i)
    
    // Update Position
    x_i = x_i + v_i
   }
   else
   {
    x_i = chaosmap(x_i)
   }
    
    // Evaluation
    If f(x_i) < f(p_best_i):
      p_best_i = x_i
    If f(p_best_i) < f(g_best):
      g_best = p_best_i
END
```

## Implementation of the code


The implementation is designed with a modular architecture that separates the functions and the maps from the particles.

### 1.  `Point`

The `Point` class represents a vector in the dominion of the problem, and it's used to represent the position of the particles.
This class exploits C++20 concepts (`Addable`, `Multipliable`, `TriComparable`) to ensure type safety for template arguments.

* **Vector Algebra:** It supports element-wise arithmetic operations (`+`, `-`, `*`, `/`) and mathematical transformations (`sin`, `cos`, `abs`, `pow`).
* **Norms:** Includes methods for calculating the 1-norm, 2-norm, and squared 2-norm.
* **Clamping:** Provides a `clamp()` method to constrain coordinates within a cube (defined by points `a` and `b`).

### 2. The Particle System

Polymorphism is exploited to manage the two different types of particles. The base `Particle` abstract class manages the `position`, `personal_best` coordinates, and the `personal_best_value`.

#### `NormalParticle` (Standard PSO)
This class implements standard Particle Swarm Optimization physics.
* **Inertia:** Uses a dynamic inertia weight $w$ that decreases linearly from 0.9 to 0.4 over the course of the iterations.
* **Coefficients:** Acceleration coefficients $c_1$ and $c_2$ are fixed at 2.5.
* **Update Logic:** The velocity is updated using the CHOPSO formula:
  $v(t+1) = w \cdot v(t) + k_1(p_{best} - x(t)) + k_2(g_{best} - x(t))$

#### `ChaoticParticle`
This class implement the chaotic version of the particle.
* **Update of the position:** Instead of calculating velocity, it updates its position by passing its current coordinates through a `ChaosMap`.

### 3. Chaos Mapping: `ChaosMap`

The `ChaosMap` class handles the mathematical transformations required for `ChaoticParticle` instances.

* **Domain Transformation:** It automatically maps coordinates between the **Global Domain** (the search space boundaries) and the **Local Domain** (the mathematical range required by the specific chaotic map).
* **Generator:** It accepts a `MapFunction` (std::function) to define the specific chaotic equation used for position generation.

### 4. Controller: `Swarm`

The `Swarm` class acts as the controller class.

* **Particles:** Keep track of all the particles (chaotic or normal) using a vector.
* **Global attributes:** Tracks the `global_best` position and value across the entire population.
* **Execution Cycle:** The `updateEveryone()` method:
    1. Triggers position updates for all particles.
    2. Updates personal bests.
    3. Recalculates the global best.

### 5. Objective Functions

It's defined a generic `Function` interface.
The actual implemented benchmark functions, which inherits from the `Function` interface are the following:

* `SphereFunction`
* `EllipsoidFunction`
* `QuinticFunction`
* `DropwaveFunction`
* `Alpine1Function`
* `AckleyFunction`



