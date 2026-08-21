# Evidence 2.8

We provide evidence for [main_paper; Conjecture 2.4](@cite) which concerns $T$-spaces of the form $X\times\mathbb{C}^2$ where $X$ is a 3d [strip geometry](@ref "Strip geometries").

## $\mathcal{A}_3\times \mathbb{C}^3$

The following Julia code verifies [main_paper; Conjecture 2.4](@cite) for $\mathcal{A}_r\times \mathbb{C}^3$ with $r=3$ up to genus 3.

The curve classes `b,c,d` correspond to the three $T$-stable rational curves  denoted and are denoted $\beta_1,\beta_2,\beta_3$ in the paper.
Note that [`get_Omega_beta`](@ref) checks all curve classes passed to it in the second argument, as well as all curve classes whose multiples appear in that list.
For example, the command `get_Omega_beta(G, [6*b], max_genus; check_predictions=true)` checks the conjecture for `b, 2*b, 3*b, 6*b` (see [Pipeline](@ref) for more details).

Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as a comment.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "Ar_C3.jl"), String) * "\n```")
```

## $X_3^{(-1,-1)}\times\mathbb{C}^2$

The following Julia code checks the conjecture up to genus 2 for $X\times\mathbb{C}^2$ where $X=X_3^{(-1,-1)}$ is the 3d strip geometry given by [`minus_one_minus_one_chain`](@ref), that is, a chain of $3$ local $\mathbb{P}^1$'s.

Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as a comment.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "minus_1_minus_1_chain.jl"), String) * "\n```")
```

## A mixed example

Finally, the following Julia code verifies [main_paper; Conjecture 2.4](@cite)
up to genus 3 for the mixed example introduced in [main_paper; Evidence 2.8](@cite).
Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as a comment.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "mixed_3d_strip.jl"), String) * "\n```")
```