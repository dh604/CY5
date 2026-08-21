# Evidence 2.13

We provide code to reproduce our evidence on [main_paper; Conjecture 2.11](@cite) which concerns $T$-spaces of the form $X\times\mathbb{C}^2$ where $X$ is the 3d [closed vertex](@ref spaces_closed_vertex).

The following Julia code verifies the conjecture [main_paper; Conjecture 2.11](@cite) up to genus 3 in the classes stated in the paper.
In the code, the curve classes `b1,b2,b3` correspond to $Q_1,Q_2,Q_3$ in the paper.

Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as a comment.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "cl_vertex.jl"), String) * "\n```")
```
