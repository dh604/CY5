# Pipeline

On this page, we explain how to analyze a given CY5 space.

## Mathematical setup

We recall the mathematical setup of [main_paper](@cite), to which we refer for full details.
Let $X$ be a smooth equivariantly Calabi-Yau GKM space with respect to the torus $T$, and let $\epsilon_1,\dots,\epsilon_r$ be the equivariant parameters for $T$.
Given a curve class $\beta$ and genus $g\ge 0$, let

```math
GW_{g,\beta}(X,T) \in \mathbb{Q}(\epsilon_1,\dots,\epsilon_r)
```

be the associated $T$-equivariant Gromov-Witten invariant (with connected domain curves) for $X$ in genus $g$ and curve class $\beta$.
Let $u$ and $Q$ be formal variables.
Then

```math
\Omega_\beta(\epsilon, u)\in \mathbb{Q}(\epsilon_1,\dots,\epsilon_r)[\![u]\!]
```

is defined by

```math
\sum_{\beta\neq 0}\sum_{g\ge 0} Q^\beta u^{2g-2}  GW_{g,\beta}(X,T)
=
\sum_{\beta\neq 0}
\sum_{k>0} \frac{1}{k} Q^{k\beta} \Omega_\beta(\epsilon,ku)
```

The point of this package is to test conjectural formulae for $\Omega_\beta$ in specific examples.

!!! note
    The conjectural formulae for $\Omega_\beta$ are formulated in terms of the substitution $q_i=e^{u\epsilon_i}$.
    As [GKMtools.jl](@cite HM25_GKMtools) computes $GW_{g,\beta}(X,T)$ using the equivariant parameters $\epsilon_i$, we only perform this substitution at the end when we compare the conjectural formula to the computed low-order terms of $\Omega_\beta$ in $u$.

## Implementation

Our pipeline has the following structure.

1. Define the $T$-space $X$ using [GKMtools.jl](@cite HM25_GKMtools).
2. Define the conjectured values of $\Omega_\beta$.
3. Compute $GW_{g,\beta}(X,T)$ using [GKMtools.jl](@cite HM25_GKMtools) for small values of $g,\beta$.
4. Solve the above recursion to compute $\Omega_\beta$ up to some order of $u$.
5. Compare the results.

Steps 1. and 2. are implemented as one function for each family of spaces that we check in the paper (see [Spaces](spaces.md)).
Steps 3.-5. are implemented in the function [`get_Omega_beta`](@ref) below.

### The main function

The following function is the main function used in our numerical experiments to confirm conjectural values of $\Omega_\beta$ up to a certain genus.

```@docs
get_Omega_beta
```

### Most important helper functions

The two most important helper functions used by [`get_Omega_beta`](@ref) are:

- [`get_GW_beta`](@ref): uses [GKMtools.jl](https://mgemath.github.io/GKMtools.jl/) to compute individual Gromov-Witten invariants and packages them in a truncated generating series.
- [`get_Omega_prediction`](@ref): provides the conjectured values of $\Omega_\beta$ for comparison to the actual values computed from Gromov-Witten invariants.

```@docs
get_GW_beta
get_Omega_prediction
```

### More details
A more detailed example of the pipeline in action can be found on the page [Details of the Pipeline](example.md).
