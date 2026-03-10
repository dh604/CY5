# Overview of Conjectures

In this section, we present for each conjecture and family of examples in [main_paper](@cite) the list of cases that we verified.

!!! note
    TODO: improve/elaborate/fill in `P3_m2_m2.jl` (as Section 7.3), `mixed_3d_strip.jl` (in Conjecture 4.2) and `free_P1_chain.jl` (as Section 4.3).
    Also fill in Conjecture 5.1 and expand material in Sections 7 and 8.

Currently, all cases are checked for some small curve classes in genus $g\le 2$.

| Section  | Space | Main function |
| ----------- | ----------- | --- |
| [Conjecture 4.2](@ref) | $X\times\mathbb{C}^2$ where $X$ is a 3d strip geometry  | [`gkm_5d_strip`](@ref) |
| [Section 4.3](@ref) | $X\times\mathbb{C}^2$ where $X$ is not a 3d strip geometry  | [`gkm_5d_p1_chain`](@ref) |
| [Conjecture 5.1](@ref) | $X_{N,m}\times\mathbb{C}^2$ (GW-gauge correspondence)  | _TODO_|
| [Conjecture 6.1](@ref) | $X\times\mathcal{A}_r$ where $X$ is a CY 3-fold  | [`X_times_Ar`](@ref) |
| [Conjecture 7.1](@ref) | $X\times\mathbb{C}^2$ where $X$ is the closed vertex  | [`gkm_5d_closed_vertex`](@ref) |
| [Section 7.2](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^2}(-1)^{\oplus 3})$  | [GKMtools.jl](@cite HM25_GKMtools) |
| [Section 7.3](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-2)^{\oplus 2})$  | [GKMtools.jl](@cite HM25_GKMtools) |
| [Conjecture 8.1](@ref) | $\mathbb{C}\times\left(\text{Tot}\,\mathcal{O}_S(-D_1)\oplus \mathcal{O}_S(-D_1)\right)$  | [`CY5_from_CY4`](@ref) |