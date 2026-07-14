
@doc raw"""
    get_Omega_beta(G::AbstractGKM_graph, betas::Vector, gMax::Int64; check_predictions::Bool=false)

Return a dictionary assigning to each curve class $\beta$ in `betas` the invariant
$\Omega_b$ as rational function in the equivariant parameters and the genus parameter $u$.

## Arguments
* `G`: the GKM graph of the space $X$ whose $\Omega_\beta$ should be computed.
* `betas`: the vector of curve classes on $X$ for which $\Omega_\beta$ should be computed.
* `gMax`: the maximum genus up to which each $\Omega_\beta$ will be approximated. That is,
          the highest order term will be $u^{2\text{gMax}-2}$.
* `check_predictions`: this optional argument with default value `false` controls whether the computed values of
          $\Omega_\beta$ are compared to the conjectured values. If `check_predictions=true`, the conjectured values
          are retrieved from [`get_Omega_prediction`](@ref).
* `show_bar`: this optional argument can be set to `false` to disable the progress bar of the Gromov-Witten localization computation.

## Example
Let us compute the genus zero and genus one terms of
$\Omega_\beta$ for the example $Z=\mathcal{A}_2\times \mathbb{C}^2$ (see [Spaces](spaces.md)).


```jldoctest pipeline_example
julia> G = Ar_times_C3(2);

julia> b = curve_class(G, Edge(1, 2))
(1, 0)

julia> c = curve_class(G, Edge(2, 3))
(0, 1)

julia> max_genus = 1;

julia> get_Omega_beta(G, [3*b, 4*b], max_genus; show_bar=false)
Calculating b=(3, 0), g=0
Calculating b=(3, 0), g=1
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(4, 0), g=0
Calculating b=(4, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 4 entries:
  (4, 0) => 0
  (2, 0) => 0
  (3, 0) => 0
  (1, 0) => (1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```
The equivariant parameters $\epsilon_1,\epsilon_2,\dots$ are displayed as `t1`, `t2`, etc.

Setting `check_predictios=true`, we can automatically check if $\Omega_\beta$ matches the output of
[`get_Omega_prediction`](@ref):

```jldoctest pipeline_example
julia> get_Omega_beta(G, [3*b, 4*b], max_genus; show_bar=false, check_predictions=true)
Calculating b=(3, 0), g=0
Calculating b=(3, 0), g=1
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(4, 0), g=0
Calculating b=(4, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1
Prediction holds for (4, 0)
Prediction holds for (2, 0)
Prediction holds for (3, 0)
Prediction holds for (1, 0)
All predictions hold.
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 4 entries:
  (4, 0) => 0
  (2, 0) => 0
  (3, 0) => 0
  (1, 0) => (1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

```
"""
function get_Omega_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; check_predictions::Bool=false, show_bar::Bool=true) where T <: CC
  
  res = get_GW_beta(G, betas, gMax; show_bar=show_bar)
  res = cc_mobius(res)

  prediction_tests = Dict{CC, Bool}()

  if check_predictions

    has_CY_substitution = has_attribute(G, :equiCY_substitution)

    if has_CY_substitution
      CY_subst = vcat([evaluate(x, t) for x in get_attribute(G, :equiCY_substitution)], [u])
    end
    
    for b in keys(res)
      p = get_Omega_prediction(G, t, u, b, gMax, res[b])

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
Return the leading part
```math
  \sum_{g \ge 0}^{\text{gMax}}u^{2g-2} GW_{g,\beta}(Z,T)
```
of $GW_{\beta}(Z,T)$ for each curve class $\beta$ listed in `beta`.

## Arguments
* `G`: the GKM graph of the GKM space $(Z,T)$ whose $GW_\beta(Z,T)$ should be computed.
* `betas`: the vector of curve classes on $Z$ for which $GW_\beta(Z,T)$ should be computed.
* `gMax`: the maximum genus up to which each $GW_\beta(Z,T)$ will be approximated. That is,
          the highest order term will be $u^{2\text{gMax}-2}$.
* `show_bar`: this optional argument can be set to `false` to disable the progress bar of the Gromov-Witten localization computation.


# Example

Let us compute $GW_\beta(Z,T)$ in degrees $\beta=1,2,3$ for $Z=\mathcal{A}_1\times\mathbb{C}^3$ as obtained using [`Ar_times_C3`](@ref).
We compute the answers up to genus 2. That is, the highest power of $u$ is $u^2$.

```jldoctest
julia> G = Ar_times_C3(1);

julia> beta = curve_class(G, Edge(1, 2));

julia> gMax = 1;

julia> get_GW_beta(G, [b, 2*b, 3*b], gMax; show_bar=false)
Calculating b=(1), g=0
Calculating b=(1), g=1
Calculating b=(2), g=0
Calculating b=(2), g=1
Calculating b=(3), g=0
Calculating b=(3), g=1
Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 3 entries:
  (2) => (1//24*t2^2*t4*u^2 + 1//24*t2^2*t5*u^2 + 1//24*t2*t4^2*u^2 + 1//24*t2*t4*t5*u^2 + 1//24*t2*t5^2*u^2 - 1//8*t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
  (1) => (1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
  (3) => (1//36*t2^2*t4*u^2 + 1//36*t2^2*t5*u^2 + 1//36*t2*t4^2*u^2 + 1//36*t2*t4*t5*u^2 + 1//36*t2*t5^2*u^2 - 1//27*t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```
"""
function get_GW_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; show_bar::Bool=true) where T <: CC
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

  return res
end