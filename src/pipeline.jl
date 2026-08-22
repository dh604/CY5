
@doc raw"""
    get_Omega_beta(G::AbstractGKM_graph, betas::Vector, gMax::Int64; check_predictions::Bool=false)

Return a dictionary assigning to each curve class $\beta$ in `betas` the invariant
$\Omega_\beta$ as rational function in the equivariant parameters and the genus parameter $u$.
Internally, this function uses [`get_GW_beta`](@ref) to get the Gromov-Witten invariants required
to compute $\Omega_\beta$.

The output also contains $\Omega_\beta$ for every $\beta$ for which some positive integer multiple is contained in `beta`.

## Arguments
* `G`: the GKM graph of the space $X$ whose $\Omega_\beta$ should be computed.
* `betas`: the vector of curve classes on $X$ for which $\Omega_\beta$ should be computed.
* `gMax`: the maximum genus up to which each $\Omega_\beta$ will be approximated. That is,
          the highest order term will be $u^{2\text{gMax}-2}$.
* `check_predictions`: this optional argument with default value `false` controls whether the computed values of
          $\Omega_\beta$ are compared to the conjectured values. If `check_predictions=true`, the conjectured values
          are retrieved from [`get_Omega_prediction`](@ref).
* `show_bar`: this optional argument can be set to `false` to disable the progress bar of the Gromov-Witten localization computation.

## Example 1: Computing $\Omega_\beta$
Let us compute the genus zero and genus one terms of
$\Omega_\beta$ for the example $Z=\mathcal{A}_2\times \mathbb{C}^3$ (see [Spaces](spaces.md)), for some examples of $\beta$.


```jldoctest pipeline_example
julia> G = Ar_times_C3(2);

julia> b = curve_class(G, Edge(1, 2))
(1, 0)

julia> c = curve_class(G, Edge(2, 3))
(0, 1)

julia> max_genus = 1;

julia> Omega = get_Omega_beta(G, [3*b, 4*b, b+c], max_genus; show_bar=false);
Calculating b=(3, 0), g=0
Calculating b=(3, 0), g=1
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(4, 0), g=0
Calculating b=(4, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1
Calculating b=(1, 1), g=0
Calculating b=(1, 1), g=1

julia> Omega[b]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

julia> Omega[2*b]
0

julia> Omega[3*b]
0

julia> Omega[4*b]
0

julia> Omega[b+c]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```
The equivariant parameters $\epsilon_1,\epsilon_2,\dots$ are displayed as `t1`, `t2`, etc.

# Example 2: Checking conjectured values of $\Omega_\beta$
Setting `check_predictios=true`, we can automatically check if $\Omega_\beta$ matches the output of
[`get_Omega_prediction`](@ref):

```jldoctest pipeline_example
julia> get_Omega_beta(G, [3*b, 4*b, b+c], max_genus; show_bar=false, check_predictions=true);
Calculating b=(3, 0), g=0
Calculating b=(3, 0), g=1
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(4, 0), g=0
Calculating b=(4, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1
Calculating b=(1, 1), g=0
Calculating b=(1, 1), g=1
Prediction holds for (4, 0)
Prediction holds for (2, 0)
Prediction holds for (3, 0)
Prediction holds for (1, 1)
Prediction holds for (1, 0)
All predictions hold.
```
In this case, all predicticted values of $\Omega_\beta$ hold in the tested range.

# Example 3: Behavior when the conjecture fails

Let us also see an example in which the conjectured value of $\Omega_\beta$ does not agree with the true value computed from Gromov-Witten invariants.
We take
```math
Z=\text{Tot}_{\mathbb{P}^1\times\mathbb{P}^1}(\mathcal{O}(-2,0)\oplus\mathcal{O}(0,-2)\oplus\mathcal{O})
```
and let $\beta$ be the curve class of one of the two factor $\mathbb{P}^1$s.
Then the properness assumption of [main_paper; Conjecture 5.1](@cite) does not hold for $\beta$ and its multiples.
Let us denote $\beta$ by `b` in the following code.

```jldoctest
julia> P1_P1 = hirzebruch_surface(NormalToricVariety, 0)
Normal toric variety

julia> P = picard_group(P1_P1)
Z^2

julia> L1 = toric_line_bundle(P1_P1, -2*P[1]) # O(-2, 0)
Toric line bundle on a normal toric variety

julia> L2 = toric_line_bundle(P1_P1, -2*P[2]) # O(0, -2)
Toric line bundle on a normal toric variety

julia> CY4 = total_space(gkm_vector_bundle_of_toric([L1, L2]))
GKM graph with 4 nodes, valency 4 and axial function:
2 -> 1 => (1, 0, 0, 0, -1, 0)
3 -> 2 => (0, 1, 0, 0, 0, -1)
4 -> 1 => (0, 1, 0, 0, 0, -1)
4 -> 3 => (-1, 0, 0, 0, 1, 0)
Standalone flags:
1.3 => (0, 0, -1, 0, -2, 0)
1.4 => (0, 0, 0, -1, 0, -2)
2.3 => (-2, 0, -1, 0, 0, 0)
2.4 => (0, 0, 0, -1, 0, -2)
3.3 => (-2, 0, -1, 0, 0, 0)
3.4 => (0, -2, 0, -1, 0, 0)
4.3 => (0, 0, -1, 0, -2, 0)
4.4 => (0, -2, 0, -1, 0, 0)

julia> G = CY5_from_CY4(CY4; equiCY=true)
GKM graph with 4 nodes, valency 5 and axial function:
2,1 -> 1,1 => (1, 0, 0, 0, -1, 0, 0)
3,1 -> 2,1 => (0, 1, 0, 0, 0, -1, 0)
4,1 -> 1,1 => (0, 1, 0, 0, 0, -1, 0)
4,1 -> 3,1 => (-1, 0, 0, 0, 1, 0, 0)
Standalone flags:
1,1.3 => (0, 0, -1, 0, -2, 0, 0)
1,1.4 => (0, 0, 0, -1, 0, -2, 0)
1,1.5 => (1, 1, 1, 1, 1, 1, 0)
2,1.3 => (-2, 0, -1, 0, 0, 0, 0)
2,1.4 => (0, 0, 0, -1, 0, -2, 0)
2,1.5 => (1, 1, 1, 1, 1, 1, 0)
3,1.3 => (-2, 0, -1, 0, 0, 0, 0)
3,1.4 => (0, -2, 0, -1, 0, 0, 0)
3,1.5 => (1, 1, 1, 1, 1, 1, 0)
4,1.3 => (0, 0, -1, 0, -2, 0, 0)
4,1.4 => (0, -2, 0, -1, 0, 0, 0)
4,1.5 => (1, 1, 1, 1, 1, 1, 0)

julia> b = curve_class(G, Edge(1, 2))
(0, 1)

julia> Omega = get_Omega_beta(G, [b, 2*b, 3*b], 1; check_predictions = true, show_bar=false);
Calculating b=(0, 1), g=0
Calculating b=(0, 1), g=1
Calculating b=(0, 2), g=0
Calculating b=(0, 2), g=1
Calculating b=(0, 3), g=0
Calculating b=(0, 3), g=1
Prediction holds for (0, 2)
Prediction holds for (0, 3)
Prediction fails for (0, 1):
  prediction = 0

julia> Omega[b]
(1//12*t1^2*t2*u^2 + 1//12*t1^2*t4*u^2 + 1//12*t1^2*t6*u^2 + 1//12*t1*t2^2*u^2 + 1//6*t1*t2*t3*u^2 + 1//6*t1*t2*t5*u^2 - 1//6*t1*t2*t6*u^2 + 1//6*t1*t3*t4*u^2 + 1//6*t1*t3*t6*u^2 + 1//6*t1*t4*t5*u^2 + 1//6*t1*t5*t6*u^2 + 1//12*t1*t6^2*u^2 - t1 + 1//12*t2^2*t3*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t3^2*u^2 + 1//6*t2*t3*t5*u^2 - 1//6*t2*t3*t6*u^2 + 1//12*t2*t5^2*u^2 - 1//6*t2*t5*t6*u^2 + 1//12*t3^2*t4*u^2 + 1//12*t3^2*t6*u^2 + 1//6*t3*t4*t5*u^2 + 1//6*t3*t5*t6*u^2 + 1//12*t3*t6^2*u^2 - t3 + 1//12*t4*t5^2*u^2 + 1//12*t5^2*t6*u^2 + 1//12*t5*t6^2*u^2 - t5)//(t1*t2*t4*u^2 + 2*t1*t2*t6*u^2 + 1//2*t1*t4^2*u^2 + t1*t4*t6*u^2 + t2^2*t4*u^2 + 2*t2^2*t6*u^2 + t2*t3*t4*u^2 + 2*t2*t3*t6*u^2 + 3//2*t2*t4^2*u^2 + t2*t4*t5*u^2 + 4*t2*t4*t6*u^2 + 2*t2*t5*t6*u^2 + 2*t2*t6^2*u^2 + 1//2*t3*t4^2*u^2 + t3*t4*t6*u^2 + 1//2*t4^3*u^2 + 1//2*t4^2*t5*u^2 + 3//2*t4^2*t6*u^2 + t4*t5*t6*u^2 + t4*t6^2*u^2)

julia> Omega[2*b]
0

julia> Omega[3*b]
0
```
As we see, the conjectured value for $\Omega_\beta$ is zero, but the true value is non-zero.
The printed output highlights whenever the conjectured value is incorrect.

Mathematically, this example is not a counterexample to [main_paper; Conjecture 5.1](@cite)
because one of the conjecture's assumptions is not satisfied in this case.
"""
function get_Omega_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; check_predictions::Bool=false, show_bar::Bool=true) where T <: CC
  
  res, t, u = get_GW_beta(G, betas, gMax; show_bar=show_bar, return_t_u=true)

  res = cc_mobius(res)

  prediction_tests = Dict{CC, Bool}()

  if check_predictions

    has_CY_substitution = has_attribute(G, :equiCY_substitution)

    if has_CY_substitution
      CY_subst = vcat([evaluate(x, t) for x in get_attribute(G, :equiCY_substitution)], [u])
    end
    
    for b in keys(res)
      # p = get_Omega_prediction(G, t, u, b, gMax, res[b])
      p = get_Omega_prediction(G, t, u, b, gMax)

      if has_CY_substitution
        p = evaluate(p, CY_subst)
      end

      prediction_tests[b] = p == res[b]
      if !prediction_tests[b]
        println("Prediction fails for $b:")
        println("  prediction = $p")
        # println("      actual = $(res[b])")
        get_attribute(G, :example_type) == :gkm_5d_strip_from_3d_CY && println("        type = $(get_attribute(G, :prediction_data)[b])")
      else
        println("Prediction holds for $b")
      end
    end
    if all(b -> prediction_tests[b], keys(prediction_tests))
      println("All predictions hold.")
    else
      println()
    end
  end

  return res
end

@doc raw"""
    get_GW_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; show_bar::Bool=true) where T <: CC

Let $(Z,T)$ be a GKM space with GKM graph $G$.
Return the truncated Grommov--Witten series
```math
\sum_{g \ge 0}^{\text{gMax}}u^{2g-2} GW_{g,\beta}(Z,T)
```
for each curve class $\beta$ listed in `beta`.
  The individual Gromov-Witten invariants are computed using the function [`gromov_witten`](https://mgemath.github.io/GKMtools.jl/stable/GW/GW/#GKMtools.gromov_witten)
  of [GKMtools](https://mgemath.github.io/GKMtools.jl/).

The output also contains the truncated Gromov-Witten series for each curve class with some positive integer multiple contained in `betas`.

## Arguments
* `G`: the GKM graph of the GKM space $(Z,T)$ whose $GW_\beta(Z,T)$ should be computed.
* `betas`: the vector of curve classes on $Z$ for which $GW_\beta(Z,T)$ should be computed.
* `gMax`: the maximum genus up to which each $GW_\beta(Z,T)$ will be approximated. That is,
          the highest order term will be $u^{2\text{gMax}-2}$.
* `show_bar`: this optional argument can be set to `false` to disable the progress bar of the Gromov-Witten localization computation.


# Example

Let us compute the Gromov-Witten series in degrees $\beta=1,2,3$ for the example

```math
Z = \text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times\mathbb{C}^3
```
and the $T$-action defined in [main_paper; Example 1.1](@cite).
Its GKM graph can be obtained using [`gkm_graph_of_example_1_1`](@ref).

We compute the Gromov-Witten series up to genus 1.
That is, the powers of $u$ that appear are $u^{-2}$ and $1$.

```jldoctest
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
"""
function get_GW_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; show_bar::Bool=true, return_t_u::Bool=false) where T <: CC
  all_betas = downward_close_ccs(betas)
  res = Dict{CC, Any}()

  # Add u variable for mobius transformation!
  S, t, u = polynomial_ring(QQ, ["t$i" for i in 1:rank_torus(G)], ["u"])
  u = u[1]

  for b in all_betas
    println("Calculating b=$b, g=0")
    gw0 = gromov_witten(G, b, 0, class_one(); g=0, show_bar=show_bar) // 1
    tmp = evaluate(gw0, t)//(u^2)
    for g in 1:gMax
      println("Calculating b=$b, g=$g")
      gw = gromov_witten(G, b, 0, class_one(); g=g, show_bar=show_bar) // 1
      tmp += evaluate(gw, t) * u^(2*g - 2)
    end
    res[b] = tmp
  end

  return return_t_u ? (res, t, u) : res
end