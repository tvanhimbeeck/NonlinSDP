## NonlinSDP
Interior-point solver for convex nonlinear optimisation over the semidefinite cone.

Copyright (C) 2024 Thomas Van Himbeeck (Licence: GLPv3)

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
where $`X\in \mathrm{H(d)}`$ is a hermitian matrix of dimension $`d`$ of subject to *semidefinite constraints* and $` f(X)`$ is a convex *concordant* real matrix function defined on the positive semidefinite cone.

### Function library
Library of *concordant* matrix functions 
| Function | formula | concavity | condition |
| -------- |-------- | --------- | --------- |
| von Neumann entropy | $`S(X) = \mathrm{tr}[ X \log(X)]`$  | concave | |
| trace function | $`t(X) = \mathrm{tr}[ f(X)]`$ | convex | |
| keyrate function    | $`h(X) = S(X) - \sum_{i} S(K_{i} X K_{i}^\dagger)`$ | convex| $`\sum_{i} K_{i}^\dagger K_{i} = \mathbf{1}`$|
| keyrate Renyi entropy | $`q_{\alpha}(X) = \sum_{i} \mathrm{tr}[(K_{i} X^{\frac{1}{\alpha}} K_{i}^\dagger)^\alpha] `$|concave| $`\sum_{i} K_{i}^\dagger K_{i} = \mathbf{1}`$ |


    A real convex or concave matrix functions is *concordant* if it satisfies the following third order condition
    ```math
    d^3 f(X)[V] \leq M |V|_X \cdot |d^2 f(X)[V]|
    ```
    for some known constant $`M`$ and for all $`X \in P(d)\succ 0, V\in H(d)`$, where $`d^k f(X)[V]`$ is the $k$th directional (Frechet) derivative and $`|V|_X = ||X^{-\frac{1}{2}} V X^{-\frac{1}{2}}||_2`$. The set of convex (or concave) *concordant* functions is closed under addition, muliplication by a positive constant, and the transformation $f(X)\mapsto f(AXA^\dagger)$.


### Solvers 

| Type | Convergence | requirements | reference |
| -- | -- | -- | -- |
|Interior-point | super-exponential $o(\log(1/\epsilon))$ | first and second order derivative, concordance property | [1,2] |
|Frank-Wolve | polynomial $O(polylog(1/\epsilon))$ | first  derivative and <a href="http://cvxr.com/cvx/">CVX <a> package |  |

### Code sample
(in progress)

### Benchmarking
(in progress)

### History
June 2022: start development
January 2024: public release

### References
1. T. Van Himbeeck (in preparation)
1. Y. Nesterov, Lectures on convex optimisation
