# NonlinSDP
Interior-point solver for convex nonlinear optimisation over the semidefinite cone.

Copyright (C) 2024 Thomas Van Himbeeck (Licence: GLPv3)

## Scope

### Convex nonlinear semidefinite programming 
The present algorithm finds solutions to convex optimisation problems of the form
```math
\begin{align}
            \text{minimize}_{X} : &\  f(X)\\
            \text{subject to} :   &\ \mathrm{tr}[X A_i] = b_i\\
                                  &\ \mathrm{tr}[X C_j] \leq d_j\\
                                  &\ X \succeq 0
\end{align}
```
where 
- $`X\in \mathrm{H(d)}`$ is a hermitian matrix of dimension $`d`$ of subject to *semidefinite constraints*
- $` f(X):\mathrm{P(d)} \rightarrow \mathbb{R}`$ is a convex *concordant* real matrix function defined on the positive semidefinite cone

### Concordant function library
This library provides the following a library of compatible functions

| Function | formula | concavity | condition |
| -------- |-------- | --------- | --------- |
| von Neumann entropy | $`S(X) = \mathrm{tr}[ X log(X)]`$  | concave | |
| trace function | $`t(X) = \mathrm{tr}[ f(X)]`$ | convex | $`f(x)`$ is convex|
| keyrate function    | $`h(X) = S(X) - \sum_{i} S(K{i} X K_{i}^\dagger)`$ | convex| $`\sum_{i} K_{i}^\dagger K_{i} = \mathbf{1}`$|
| keyrate Renyi entropy | $`q_{\alpha}(X) = \sum_{i} \mathrm{tr}[(K_{i} X^{\frac{1}{\alpha}} K_{i}^\dagger)^\alpha] `$|concave| $`\sum_{i} K_{i}^\dagger K_{i} = \mathbf{1}`$ |

Any concave/convex matrix function satisfy the following *concordance* property can be added:
```math
d^3 f(X)[V] \leq M |V|_X \cdot |d^2 f(X)[V]|
```
for some known constant $`M`$ and for all $`X \in P(d)\succ 0, V\in H(d)`$, where $`d^k f(X)[V]`$ is the $k$th directional (Frechet) derivative and $`|V|_X = ||X^{-\frac{1}{2}} V X^{-\frac{1}{2}}||_2`$.

The concordance property is closed under addition, muliplication by a positive constant, and the transformation $f(X)\mapsto f(AXA^\dagger)$.

## Solvers 

| Type | Convergence | requirements | reference |
| -- | -- | -- | -- |
|(Path-following) interior-point | super-exponential $o(\log(1/\epsilon))$ | first and second order derivative, concordance property | [1,2] |
|Frank-Wolve | polynomial $O(polylog(1/\epsilon))$ | first  derivative and <a href="http://cvxr.com/cvx/">CVX <a> package |  |

### Benchmarking

### References
1. T. Van Himbeeck (in preparation)
1. Y. Nesterov, Lectures on convex optimisation
