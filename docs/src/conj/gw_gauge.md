# Evidence 3.5

On this page, we provide Julia code to verify [main_paper; Conjecture 3.3](@cite) concerning the $T$-space $X_{N,m}\times\mathbb{C}^2$ obtained from [`gkm_5d_gauge`](@ref).
The curve classes `F[1],...,F[N-1]` and `B[1],...,B[N]` in the code correspond to $F_1,\dots,F_{N-1}$ and $C_1,\dots,C_N$ in the paper, respectively, see [main_paper; Figure 7](@cite).

Under each line calling [`get_Omega_beta`](@ref), we recorded the crucial line of the output as a comment.

## $N=2, m=-1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_m1.jl"), String) * "\n```")
```

## $N=2, m=0$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_0.jl"), String) * "\n```")
```

## $N=2, m=1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N2_1.jl"), String) * "\n```")
```

## $N=3, m=-2$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m2.jl"), String) * "\n```")
```

## $N=3, m=-1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m1.jl"), String) * "\n```")
```

## $N=3, m=0$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_0.jl"), String) * "\n```")
```

## $N=3, m=1$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_1.jl"), String) * "\n```")
```

## $N=3, m=2$

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_2.jl"), String) * "\n```")
```

## Necessity of $|m| < N$

Let $N=3$ and $m=4$ so that $|m|\ge N$.
In this example, all equivariant parameters `t1`, `t3`, `t4`, `t5` appear in $\Omega_\beta$, which shows that $|m|<N$ is in general a necessary assumption for [main_paper; Conjecture 3.3](@cite) to hold.
Indeed, for $|m|<N$, the conjectured formula only contains the equivariant parameters `t4` and `t5`.
(Note that `t2` does not appear because it is expressed in terms of the other parameters on the equivariantly Calabi-Yau torus.)

Under each call to [`get_Omega_beta`](@ref), we include the output of the computed values for $\Omega_\beta$ as a comment, showing the equivariant parameters `t1`, `t3`, `t4`, `t5`.

```@eval
using Markdown
Markdown.parse("```julia\n" * read(joinpath(@__DIR__, "..", "..", "..", "test", "results", "gauge", "gauge_N3_m4.jl"), String) * "\n```")
```