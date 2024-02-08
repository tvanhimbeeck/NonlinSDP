# NonlinSDP
Solver solutions for convex nonlinear semidefinite programming in MATLAB with applications to key rate calculations in  *Quantum Key Distribution* and other entropic optimisation problems.

Copyright (C) 2024 Thomas Van Himbeeck (Licence: GLPv3)

## Convex nonlinear semidefinite programming 
The present algorithm finds solutions to convex optimisation problems of the form
```math
\begin{align}
            \text{minimize}_{X}& f(X)\\
            \text{subject to}  & tr[X A_i] = b_i\\
                               & tr[X C_j] <= d_j\\
                               & X >=0
\end{align}
```
where 
- $`X\in \mathrm{H}`$ is a hermitian matrix of subject to **semidefinite constraints**
- $`f(X)`:\mathrm{H} \rightarrow \mathbb{R}$ is a convex real matrix function **compatible with the semidefinite barier**: 
```math
d^3 f(X)[V] \leq d^2 f(X)[V] |V|_X
```
for all $`X \succ 0`$

## Function library
The solver provides a library of functions

| Function | formula | concavity |
| -------- |-------- | --------- |
| von Neumann entropy | S(X) = tr[ X log(X)]  | concave |
| keyrate function    | \sum_p S(K_p X K_p) - \sum_{kp} S

Other convex functions can be added to the function library provided they satisfy a technical assumption called concordance (see below), and the derivative and hessian are given.

### Jointly-convex optimisation problems
The code also supports extentions of this problem such as jointly-convex matrix functions

            minimize_{X}   g(X,Y)
              subject to   tr[X A_i] + tr[Y B_i] = c_i
                           tr[X D_j] + tr[Y E_j] <= f_j
                           X,Y >=0
where *g(X,Y)* is a **joint-convex matrix function** and *X* and *Y* semidefinite hermiatian matrices that satisfy joint linear equality and inequality constaints.

## Interior-point solver

### Features

#### Dual-problem / Lower-bounds
Upper- and lower- bound are provided on the optimal solution, that agree up to the desired precision.

#### Similar to Semidefinite programming
The code is based on standard state-of-the-art path-following interior-point algorithm. The same type algorithms are behind the efficient solvers for Semidefinite programming.

#### Exponential convergence
The path-following interior-point algorithm behind this code guaranties exponential convergence: The running time of the code is O(log(1/epsilon)), where epsilon is the desired precision. Moreover, the expected running time can be estimated at the start as the algorithm returns an upper-bound on the number of search iteration needed to reach a given precision.

### Workings
The code is based on standard state-of-the-art path-following interior-point algorithm. A self-concordant barrier function is constructed for the epigraph of the matrix function.

The code is based on the interior-point path following.

### Improvements

It is expected that running time improvements could be obtained in the future.

### Provide code it is based on
