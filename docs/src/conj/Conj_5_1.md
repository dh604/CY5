# Conjecture 5.1

!!! note
    TODO: fill in material on Conjecture 5.1

## Tested cases

### $N=3, m=2$

In the following example with $|m|<N$, only the equivariant parameters `t4` and `t5` appear.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "gauge", "gauge_N3_m2.jl"), String) * "\n```")
```

### $N=3, m=4$

In the following example with $|m|\ge N$, all equivariant parameters `t1`, `t3`, `t4`, `t5` appear.
(Recall that on the equivariantly Calabi-Yau torus, `t2` is expressed in terms of the other parameters.)

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "gauge", "gauge_N3_m4.jl"), String) * "\n```")
```