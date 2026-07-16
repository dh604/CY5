# Julia code in the paper

On this page, we mirror the Julia code presented in the paper.

## Example 1.1

In [main_paper; Example 1.1](@cite), the GKM graph $G$ of
$Z = \text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times \mathbb{C}^3$
with a certain torus action is introduced.
The following function can be used to obtain its GKM graph.

```@docs
gkm_graph_of_example_1_1
```

## Example 1.2

In [main_paper; Example 1.2](@cite), we compute the Gromov-Witten of $Z$ in degrees $1,2,3$ and genus $0$ and $1$ using the following code.

```jldoctest Example_from_paper
julia> using Oscar, GKMtools, CY5

julia> G = gkm_graph_of_example_1_1();

julia> beta = curve_class(G, Edge(1, 2));

julia> gMax = 1;

julia> GW = get_GW_beta(G, [beta, 2*beta, 3*beta], gMax; show_bar=false);
Calculating b=(1), g=0
Calculating b=(1), g=1
Calculating b=(2), g=0
Calculating b=(2), g=1
Calculating b=(3), g=0
Calculating b=(3), g=1

julia> GW[beta]
(1//12*t1*t4*t5*u^2 + 1//12*t1*t4*t6*u^2 + 1//12*t1*t5*t6*u^2 + t1 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t4*t6*u^2 + 1//12*t2*t5*t6*u^2 + t2 + 1//12*t3*t4*t5*u^2 + 1//12*t3*t4*t6*u^2 + 1//12*t3*t5*t6*u^2 + t3)//(t4*t5*t6*u^2)

julia> GW[2*beta]
(1//24*t1*t4*t5*u^2 + 1//24*t1*t4*t6*u^2 + 1//24*t1*t5*t6*u^2 + 1//8*t1 + 1//24*t2*t4*t5*u^2 + 1//24*t2*t4*t6*u^2 + 1//24*t2*t5*t6*u^2 + 1//8*t2 + 1//24*t3*t4*t5*u^2 + 1//24*t3*t4*t6*u^2 + 1//24*t3*t5*t6*u^2 + 1//8*t3)//(t4*t5*t6*u^2)

julia> GW[3*beta]
(1//36*t1*t4*t5*u^2 + 1//36*t1*t4*t6*u^2 + 1//36*t1*t5*t6*u^2 + 1//27*t1 + 1//36*t2*t4*t5*u^2 + 1//36*t2*t4*t6*u^2 + 1//36*t2*t5*t6*u^2 + 1//27*t2 + 1//36*t3*t4*t5*u^2 + 1//36*t3*t4*t6*u^2 + 1//36*t3*t5*t6*u^2 + 1//27*t3)//(t4*t5*t6*u^2)
```

## Example 1.4

Finally, in [main_paper; Example 1.4](@cite) we compute $\Omega_\beta$ for the same example and in the same range of degree and genus.
We do this by continuing the above code as follows.

```jldoctest Example_from_paper
julia> Omega = get_Omega_beta(G, [beta, 2*beta, 3*beta], gMax; show_bar=false);
Calculating b=(1), g=0
Calculating b=(1), g=1
Calculating b=(2), g=0
Calculating b=(2), g=1
Calculating b=(3), g=0
Calculating b=(3), g=1

julia> Omega[beta]
(1//12*t1*t4*t5*u^2 + 1//12*t1*t4*t6*u^2 + 1//12*t1*t5*t6*u^2 + t1 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t4*t6*u^2 + 1//12*t2*t5*t6*u^2 + t2 + 1//12*t3*t4*t5*u^2 + 1//12*t3*t4*t6*u^2 + 1//12*t3*t5*t6*u^2 + t3)//(t4*t5*t6*u^2)

julia> Omega[2*beta]
0

julia> Omega[3*beta]
0
```