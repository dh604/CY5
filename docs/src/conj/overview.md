# Overview

In this section, we provide code to reproduce the results cited as *Evidence* in [main_paper](@cite).

!!! warning
    **TODO:** Restructure this section when the paper is ready enough.

| Section  | Space | Main function |
| ----------- | ----------- | --- |
| [Conjecture 4.2](@ref) | $X\times\mathbb{C}^2$ where $X$ is a 3d strip geometry  | [`gkm_5d_strip`](@ref) |
| [Section 4.3](@ref) | $X\times\mathbb{C}^2$ where $X$ is not a 3d strip geometry  | [`gkm_5d_p1_chain`](@ref) |
| [Conjecture 5.1](@ref) | $X_{N,m}\times\mathbb{C}^2$ (GW-gauge correspondence)  | [`gkm_5d_gauge`](@ref) |
| [Conjecture 6.1](@ref) | $X\times\mathcal{A}_r$ where $X$ is a CY 3-fold  | [`X_times_Ar`](@ref) |
| [Conjecture 7.1](@ref) | $X\times\mathbb{C}^2$ where $X$ is the closed vertex  | [`gkm_5d_closed_vertex`](@ref) |
| [Section 7.2](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^2}(-1)^{\oplus 3})$  | [`gkm_5d_P2_111`](@ref) |
| [Section 7.3](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-2)^{\oplus 2})$  | [`gkm_5d_P3_22`](@ref) |
| [Section 7.4](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-1)\oplus\mathcal{O}_{\mathbb{P}^3}(-3))$  | [`gkm_5d_P3_13`](@ref) |
| [Section 7.5](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(0)\oplus\mathcal{O}_{\mathbb{P}^3}(-4))$  | [`gkm_5d_P3_04`](@ref) |
| [Conjecture 8.1](@ref) | $\mathbb{C}\times\left(\text{Tot}\,\mathcal{O}_S(-D_1)\oplus \mathcal{O}_S(-D_1)\right)$  | [`CY5_from_CY4`](@ref) |


!!! note
    Note that we executed most of the code listed in this section on a compute cluster provided by the Mathematics Department of ETH Zürich with the following metrics.
    - **Operating System:** Fedora
    - **CPU:** 2x Intel Xeon Gold 6254 (18 Cores) 3.1 GHz
    - **RAM:** 512 GB

    Therefore, it may be hard to reproduce some of our examples on an ordinary laptop.