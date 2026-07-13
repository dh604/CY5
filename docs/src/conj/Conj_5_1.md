# Conjecture 5.1

## Tested cases

### $N=2, m=-1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_m1.jl"), String) * "\n```")
```

### $N=2, m=0$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_0.jl"), String) * "\n```")
```

### $N=2, m=1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_1.jl"), String) * "\n```")
```

### $N=3, m=-2$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m2.jl"), String) * "\n```")
```

### $N=3, m=-1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m1.jl"), String) * "\n```")
```

### $N=3, m=0$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_0.jl"), String) * "\n```")
```

### $N=3, m=1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_1.jl"), String) * "\n```")
```

### $N=3, m=2$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_2.jl"), String) * "\n```")
```

## Necessity of $|m| < N$

### $N=3, m=4$

In the following example with $|m|\ge N$, all equivariant parameters `t1`, `t3`, `t4`, `t5` appear, underlining the necessity of the assumption $|m|<N$ for our conjecture to hold.
Indeed, for $|m|<N$, the conjectured formula only contains the equivariant parameters `t4` and `t5`.
(Note that `t2` does not appear because it is expressed in terms of the other parameters on the equivariantly Calabi-Yau torus.)

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m4.jl"), String) * "\n```")
```