# Conjecture 6.1

This conjecture concerns equivariantly Calabi-Yau $T$-spaces of the form $Z=X\times\mathcal{A}_r$ where $X$ a 3-fold.

## Statement
Recall from [main_paper; Conjecture 6.1](@cite) that the conjecture states
```math
\Omega_\beta = 0
```
for any $\beta$ such that $\overline{M}_g(Z,\beta)$ is proper for all $g\ge 0$.

## Verified cases

#### The case $\text{Tot}(K_{\mathbb{P}^2})\times \mathcal{A}^r$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "KP2_A3.jl"), String) * "\n```")
```

#### The case $\text{Tot}(K_{\Sigma_5})\times \mathcal{A}^r$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "KH5_A3.jl"), String) * "\n```")
```

#### The case $\mathcal{A}_s\times\mathbb{C}\times \mathcal{A}^r$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "As_C_Ar.jl"), String) * "\n```")
```