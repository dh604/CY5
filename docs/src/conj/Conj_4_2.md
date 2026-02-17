# Conjecture 4.2

This conjecture concerns $T$-spaces of the form $X\times\mathbb{C}^2$ where $X$ is a 3d [strip geometry](@ref "Strip geometries").

**Notation:**
Let $(\epsilon_i)_i$ be the equivariant parameters for $T$, $u$ the formal variable tracking the genus, $q_i:=e^{u\epsilon_i}$, and write
```math
[q] := q^{1/2} - q^{-1/2}.
```

## Statement
Recall from [main_paper; Conjecture 4.2](@cite) that the conjecture states
```math
\Omega_\beta = -\frac{[q_2q_4q_5]}{[q_2][q_4][q_5]}
```
if $\beta$ corresponds to an element of $I_u$,
```math
\Omega_\beta = -\frac{[q_3q_4q_5]}{[q_3][q_4][q_5]}
```
if $\beta$ corresponds to an element of $I_d$, and
```math
\Omega_\beta = \frac{1}{[q_4][q_5]}
```
if $\beta$ corresponds to an element of $I_0$.

Moreover, if $\beta$ is not given by a connected chain of $T$-stable rational curves in $X$, all with multiplicity one, then
```math
\Omega_\beta=0.
```

## Verified cases

#### The case $\mathcal{A}_r\times \mathbb{C}^3$

The following code verifies the conjecture for $\mathcal{A}_r\times \mathbb{C}^3$ with $r=3$ up to genus 2.

The curve classes `b,c,d` correspond to the three $T$-stable rational curves.
Note that [`get_Omega_beta`](@ref) checks all curve classes passed to it in the second argument, as well as all curve classes whose multiples appear in that list.
For example, the command `get_Omega_beta(G, [6*b], max_genus; check_predictions=true)` checks the conjecture for `b, 2*b, 3*b, 6*b`.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "Ar_C3.jl"), String) * "\n```")
```

#### Chains of $\mathcal{O}(-1)\oplus\mathcal{O}(-1)$

The following code checks the conjecture up to genus 2 for $X\times\mathbb{C}^2$ where $X$ is a 3d strip geometry given by [`minus_one_minus_one_chain`](@ref), that is, a chain of $r$ local $\mathbb{P}^1$'s.
We set $r=3$.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "minus_1_minus_1_chain.jl"), String) * "\n```")
```