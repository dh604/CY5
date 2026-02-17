# Conjecture 7.1

This conjecture concerns $T$-spaces of the form $X\times\mathbb{C}^2$ where $X$ is the 3d [closed vertex](@ref "Closed vertex").

**Notation:**
Let $(\epsilon_i)_i$ be the equivariant parameters for $T$, $u$ the formal variable tracking the genus, $q_i:=e^{u\epsilon_i}$, and write
```math
[q] := q^{1/2} - q^{-1/2}.
```

## Statement
Let $\beta_1,\beta_2,\beta_3$ be the curve classes corresponding to $Q_1,Q_2,Q_3$ in the notation of [main_paper; Conjecture 7.1](@cite).
Recall that the conjecture states
```math
\Omega_\beta = \frac{1}{[q_4][q_5]}
```
if $\beta = \beta_i$ for some $i\in\{1,2,3\}$ or $\beta=\beta_1+\beta_2+\beta_3$, and
```math
\Omega_\beta = -\frac{[q_iq_4q_5]}{[q_i][q_4][q_5]}
```
if $\beta=\beta_1 + \beta_2 + \beta_3 - \beta_i$ for some $i\in\{1,2,3\}$.
For all other $\beta$, the conjecture states
```math
\Omega_\beta = 0
```

## Verified cases


```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "cl_vertex.jl"), String) * "\n```")
```