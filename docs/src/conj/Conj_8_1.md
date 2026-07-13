# Conjecture 8.1

This conjecture concerns equivariantly Calabi-Yau $T$-spaces of the form $Y\times\mathbb{C}$ where
```math
Y = \text{Tot}\,\mathcal{O}_S(-D_1)\oplus \mathcal{O}_S(-D_1)
```
where $S$ is a projective surface and $D_1,D_2$ divisors satisfying $D_1+D_2 \in |-K_S|$.


## Statement
The statement of [main_paper; Conjecture 8.1](@cite) is that, under a suitable positivity condition on $D_1,D_2$, we have $\Omega_\beta=0$ for all $\beta$
such that $\overline{M}_g(Y,\beta)$ is proper for all $g\ge 0$. 

## Verified cases

#### The case $\mathcal{O}(-1)\oplus \mathcal{O}(-2)$ on $\mathbb{P}^2$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P2.jl"), String) * "\n```")
```

#### The case $\mathcal{O}(-1, -1)\oplus \mathcal{O}(-1, -1)$ on $\mathbb{P}^1\times\mathbb{P}^1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P1_P1_1_1.jl"), String) * "\n```")
```


#### The case $\mathcal{O}(-2, 0)\oplus \mathcal{O}(0, -2)$ on $\mathbb{P}^1\times\mathbb{P}^1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P1_P1_0_2.jl"), String) * "\n```")
```
