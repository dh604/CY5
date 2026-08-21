# Evidence 4.8

On this page, we present Julia code to reproduce the claims made in [main_paper; Evidence 4.8](@cite) concerning the values of $\Omega_{dL}$ for
```math
\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-a)\oplus\mathcal{O}_{\mathbb{P}^3}(a-4))
```
with $a\in\{0,1,2\}$.
Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as comment.

## $\mathcal{O}_{\mathbb{P}^3}(0)\oplus \mathcal{O}_{\mathbb{P}^3}(-4)$

In this case ($a=0$), the GKM graph is obtained using [`gkm_5d_P3_04`](@ref).

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "P3_04.jl"), String) * "\n```")
```

## $\mathcal{O}_{\mathbb{P}^3}(-1)\oplus \mathcal{O}_{\mathbb{P}^3}(-3)$

In this case ($a=1$), the GKM graph is obtained using [`gkm_5d_P3_13`](@ref).

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "P3_13.jl"), String) * "\n```")
```


## $\mathcal{O}_{\mathbb{P}^3}(-2)^{\oplus 2}$

In this case ($a=2$), the GKM graph is obtained using [`gkm_5d_P3_22`](@ref).

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "P3_22.jl"), String) * "\n```")
```