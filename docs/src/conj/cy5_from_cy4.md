# Evidence 5.5

On this page, we provide Julia code to verify the vanishing of certain Gromov-Witten invariants conjectured in [main_paper; Conjecture 5.1](@cite).
The corresponding GKM graphs are obtained using [`CY5_from_CY4`](@ref).
After each line calling [`get_Omega_beta`](@ref), we added the crucial line of the output as a comment.


## $\mathcal{O}(-1)\oplus \mathcal{O}(-2)$

This example concerns the total space of $\mathcal{O}(-1)\oplus \mathcal{O}(-2)\oplus\mathcal{O}$ on $\mathbb{P}^2$.
The curve calss `b` in the code corresponds to $H$ in the paper.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P2.jl"), String) * "\n```")
```

## $\mathcal{O}(-1, -1)\oplus \mathcal{O}(-1, -1)$

This example concerns the total space of $\mathcal{O}(-1, -1)\oplus \mathcal{O}(-1, -1)\oplus\mathcal{O}$ on $\mathbb{P}^1\times\mathbb{P}^1$.
The curve classes `b,c` in the code correspond to $(1, 0)$ and $(0, 1)$ in the paper.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P1_P1_1_1.jl"), String) * "\n```")
```


## $\mathcal{O}(-2, 0)\oplus \mathcal{O}(0, -2)$


This example concerns the total space of $\mathcal{O}(-2, 0)\oplus \mathcal{O}(0, -2)\oplus\mathcal{O}$ on $\mathbb{P}^1\times\mathbb{P}^1$.
The curve classes `b,c` in the code correspond to $(1, 0)$ and $(0, 1)$ in the paper.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "CY4_vanishing", "P1_P1_0_2.jl"), String) * "\n```")
```

The final line of code does not contradict [main_paper; Conjecture 5.1](@cite).
Instead, it confirms that the properness assumption of the conjecture is necessary.