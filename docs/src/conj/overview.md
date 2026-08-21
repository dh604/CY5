# Overview of GW evidence

In this section, we explain how our GW calculations cited as *Evidence* in [main_paper](@cite) can be reproduced.

| Evidence  | Space | Main function(s) |
| ----------- | ----------- | --- |
| [Evidence 2.8](@ref) | $X\times\mathbb{C}^2$ where $X$ is a 3d strip geometry  | [`gkm_5d_strip`](@ref) |
| [Evidence 2.13](@ref) | $X\times\mathbb{C}^2$ where $X$ is the closed vertex  | [`gkm_5d_closed_vertex`](@ref) |
| [Evidence 3.5](@ref) | $X_{N,m}\times\mathbb{C}^2$ (GW-gauge correspondence)  | [`gkm_5d_gauge`](@ref) |
| [Evidence 4.2](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^2}(-1)^{\oplus 3})$  | [`gkm_5d_P2_111`](@ref) |
| [Evidence 4.8](@ref) | $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-a)\oplus\mathcal{O}_{\mathbb{P}^3}(a-4))$  | [`gkm_5d_P3_04`](@ref)/[`_13`](@ref gkm_5d_P3_13)/[`_22`](@ref gkm_5d_P3_22) |
| [Evidence 5.5](@ref) | $\mathbb{C}\times\left(\text{Tot}\,\mathcal{O}_S(-D_1)\oplus \mathcal{O}_S(-D_1)\right)$  | [`CY5_from_CY4`](@ref) |

In the last row, we treat the cases $S=\mathbb{P}^2$ with $\mathcal{O}(-1)\oplus\mathcal{O}(-2)$ and $S=\mathbb{P}^1\times \mathbb{P}^1$ with $\mathcal{O}(-1,-1)\oplus\mathcal{O}(-1,-1)$ and $\mathcal{O}(-2, 0)\oplus\mathcal{O}(0, -2)$.

!!! note
    Note that we executed most of the code listed in this section on a compute cluster provided by the Mathematics Department of ETH Zürich with the following metrics.
    - **Operating System:** Fedora
    - **CPU:** 2x Intel Xeon Gold 6254 (18 Cores) 3.1 GHz
    - **RAM:** 512 GB

    Therefore, it may be hard to reproduce some of our examples on an ordinary PC.