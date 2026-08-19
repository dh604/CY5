# GW_CY5.jl

_A Julia package to compute the membrane indices_ $\Omega_\beta$ _of toric CY5s._

We build on the Julia package [GKMtools.jl](https://mgemath.github.io/GKMtools.jl/) [HM25_paper, HM25_GKMtools](@cite) to compute equivariant Gromov-Witten invariants of the relevant spaces.

!!! warning
    **TODO:** complete read-through once the article is reasonably ready, to make all notation and references consistent with the article.

## Contents

* [Installation instructions](@ref "Installation")
* [Pipeline](pipeline.md): explanation of the workflow of analyzing any given example with this package.
* [Spaces](spaces.md): list of spaces that are implemented to work with the pipeline.
* Examples:
    - [Julia code in the paper](code_from_paper.md): a mirror of all Julia code contained in our paper.
    - [Details of the Pipeline](example.md): an example of the Pipeline, worked out in detail, showing some internal details.
* [Evidence](conj/overview.md): for each numerical experiment cited as *Evidence* in the paper, we provide code to reproduce the stated results.


## Installation

First, install **OSCAR** using the instructions [here](https://www.oscar-system.org/install/).
After that, open Julia and install **GW_CY5.jl** by typing the following code.
This will automatically install the package **GKMtools.jl** on which GW_CY5.jl builds.

```julia-repl
julia> using Pkg;
julia> Pkg.add(url="https://github.com/dh604/CY5")
```

You can now load the package by typing the following line.

```julia-repl
julia> using Oscar, GKMtools, GW_CY5
```