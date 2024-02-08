# NonlinSDP
Solver solutions for convex nonlinear semidefinite programming in MATLAB with applications to key rate calculations in  *Quantum Key Distribution* and other entropic optimisation problems.

Copyright (C) 2024 Thomas Van Himbeeck (Licence: GLPv3)

### Convex nonlinear semidefinite programming 
The present algorithm finds solutions to convex optimisation problems of the form
```math
\begin{align}
            \text{minimize}_{X} : &\  f(X)\\
            \text{subject to} :   &\ \mathrm{tr}[X A_i] = b_i\\
                                  &\ \mathrm{tr}[X C_j] <= d_j\\
                                  &\ X >=0
\end{align}
```
where 
- $`X\in \mathrm{H(d)}`$ is a hermitian matrix of dimension d of subject to **semidefinite constraints**
- $`f(X)`:\mathrm{H(d)} \rightarrow \mathbb{R}`$ is a convex real matrix **compatible with the semidefinite barier**, ie.\ such that
```math
d^2 f(X)[V] \leq |V|_X \cdot d^2 f(X)[V]  \qquad \text{, for all } 
```
for all $`X \succ 0, V\in H(d)`$ where $`d^k f(X)[V]`$ is the $`k`$th directional (Frechet) derivative and $`|V|_X = ||X^{-\frac{1}{2}} V X^{-\frac{1}{2}}||_2`$

### Function library
This library provides the following a library of functions

| Function | formula | concavity | condition |
| -------- |-------- | --------- | --------- |
| von Neumann entropy | $`S(X) = tr[ X log(X)]`$  | concave | |
| trace function | $`t(X) = \Tr[ f(X)] | convex | $`f(x)`$ is convex|
| keyrate function    | $`h(X) = S(X) - \sum_{kp} S(K{i} X K_{i}^\dagger)`$ | convex| \sum_i K_i^\dagger K_i = \mathbf 1|
| keyrate Renyi entropy | $`q_\alpha = \sum_i \mathrm{tr}[(K_i X^{\frac{1/\alpha}} K_i^\dagger)^\alpha] `$|concave|

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
