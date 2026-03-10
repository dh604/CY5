
@doc raw"""
    get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta, max_genus::Int64, res)

Return the conjectural prediction of $\Omega_\beta$ for the space with GKM graph `G`.

## Arguments
* `G`: the GKM graph of the space $X$ whose $\Omega_\beta$ we care about.
* `t`: the vector of equivariant parameters.
* `u`: the formal variable $u$ that keeps track of the genus.
* `beta`: the curve class $\beta$ for which the conjectural value of $\Omega_\beta$ should be returned.
* `max_genus`: the maximum genus up to which the conjectural value of $\Omega_\beta$ should be returned.
      That is, the highest exponent of $u$ in the returned value is `2*max_genus - 1`.

!!! warning
    The field `res` is needed as the prediction is phrased in terms of the calculated genus zero invariants
    for some cases. This will be removed in future.

!!! note
    This function requires the attribute `:example_type` of `G` to be set.
    This holds automatically whenever `G` is obtained from one of the constructors in the list
    of [implemented spaces](spaces.md).

## Example
Let us get the predicted $\Omega_\beta$ for $Z=\mathcal{A}_2\times \mathbb{C}^2$, continuing the example
$Z=\mathcal{A}_2\times \mathbb{C}^2$ from [`get_Omega_beta`](@ref).

```jldoctest pipeline_example_2
julia> G = Ar_times_C3(2);

julia> b = curve_class(G, Edge(1, 2));

julia> c = curve_class(G, Edge(2, 3));

julia> max_genus = 1;

julia> R, (t1, t2, t3, t4, t5, u) = polynomial_ring(QQ, ["t1", "t2", "t3", "t4", "t5", "u"])
(Multivariate polynomial ring in 6 variables over QQ, QQMPolyRingElem[t1, t2, t3, t4, t5, u])

julia> prediction = get_Omega_prediction(G, [t1, t2, t3, t4, t5], u, b, max_genus, nothing)
(-1//12*t3^2*t4*u^2 - 1//12*t3^2*t5*u^2 - 1//12*t3*t4^2*u^2 - 1//4*t3*t4*t5*u^2 - 1//12*t3*t5^2*u^2 - t3 - 1//12*t4^2*t5*u^2 - 1//12*t4*t5^2*u^2 - t4 - t5)//(t3*t4*t5*u^2)
```

This does not yet match the result of [`get_Omega_beta`](@ref) above, because we still need
to make the equivariantly Calabi-Yau substitution.
The appropriate substitution is stored in `G`:

```jldoctest pipeline_example_2
julia> substitution = get_attribute(G, :equiCY_substitution)
5-element Vector{QQMPolyRingElem}:
 t1
 t2
 -t2 - t4 - t5
 t4
 t5

julia> # Translate the substitution vector to our own variables.
        substitution = [evaluate(t, [t1, t2, t3, t4, t5]) for t in substitution]
5-element Vector{QQMPolyRingElem}:
 t1
 t2
 -t2 - t4 - t5
 t4
 t5

julia> evaluate(prediction, vcat(substitution, [u]))
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```
This matches the result of the example in [`get_Omega_beta`](@ref), as required.

The above is an illustration of how the pipeline works.
When [`get_Omega_beta`](@ref) is used with `check_predictions=true`, it uses the same mechanism
to determine if the actual $\Omega_\beta$ matches the conjectured one.
"""
function get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64, res)
  @req has_attribute(G, :example_type) "G has no prediction for Omega."
  et = get_attribute(G, :example_type)

  if et == :gkm_5d_strip_from_3d_CY
    return gkm_5d_strip_prediction(G, t, u, beta, max_genus)
  elseif et == :X_times_Ar
    return zero(u)
  elseif et == :closed_vertex
    return gkm_5d_closed_vertex_prediction(G, t, u, beta, max_genus)
  elseif et == :CY5_from_CY4
    return zero(u)
  elseif et == :P1_chain_5d
    return gkm_5d_free_strip_prediction(G, t, u, beta, max_genus, res)
  else
    error("Example type $et is not implemented.")
  end
end

function gkm_5d_strip_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)

  @req !iszero(beta) "beta must not be zero"

  max_deg = 2*max_genus - 2
  
  # only nonzero exponents should be one
  b2 = rank(parent(beta))
  any(i -> beta[i] > 1, 1:b2) && return zero(u)

  
  # nonzero exponents should be adjacent
  first_nonzero = findfirst(i -> !iszero(beta[i]), 1:b2)
  last_nonzero = findlast(i -> !iszero(beta[i]), 1:b2)
  any(i -> iszero(beta[i]), first_nonzero+1:last_nonzero-1) && return zero(u)

  # determine if we are in I_u, I_d, or I_0.
  strip_dict = get_attribute(G, :prediction_data)
  I_type = strip_dict[beta]

  if I_type == :I_0
    return gkm_5d_strip_prediction_I_0(t, u, max_deg)
  elseif I_type == :I_u
    return gkm_5d_strip_prediction_I_u(t, u, max_deg)
  elseif I_type == :I_d
    return gkm_5d_strip_prediction_I_d(t, u, max_deg)
  else
    error("strip_dict[$beta] has invalid value $I_type.")
  end

end

function gkm_5d_strip_prediction_I_u(t, u, max_deg::Int64)
  return -1 * t_sq_bracket_frac([u * (t[2] + t[4] + t[5])], u .* [t[2], t[4], t[5]], max_deg)
end

function gkm_5d_strip_prediction_I_d(t, u, max_deg::Int64)
  return -1 * t_sq_bracket_frac([u * (t[3] + t[4] + t[5])], u .* [t[3], t[4], t[5]], max_deg)
end

function gkm_5d_strip_prediction_I_0(t, u, max_deg::Int64)
  return t_sq_bracket_frac([], u .* [t[4], t[5]], max_deg)
end

function gkm_5d_closed_vertex_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)

  @req !iszero(beta) "beta must not be zero"

  b1 = curve_class(G, "0,1", "1,1")
  b2 = curve_class(G, "0,1", "2,1")
  b3 = curve_class(G, "0,1", "3,1")
  H2 = parent(b1)
  @req (b1 == gens(H2)[1]) && (b2 == gens(H2)[2]) && (b3 == gens(H2)[3]) "Generator labeling of H2(closed vertex) changed."

  max_deg = 2*max_genus - 2
  
  if beta in [b1, b2, b3, b1+b2+b3]
    return t_sq_bracket_frac([], u .* [t[4], t[5]], max_deg)
  elseif beta == b1+b2
    return -t_sq_bracket_frac([u * (t[3] + t[4] + t[5])], u .* [t[3], t[4], t[5]], max_deg)
  elseif beta == b2+b3
    return -t_sq_bracket_frac([u * (t[1] + t[4] + t[5])], u .* [t[1], t[4], t[5]], max_deg)
  elseif beta == b1+b3
    return -t_sq_bracket_frac([u * (t[2] + t[4] + t[5])], u .* [t[2], t[4], t[5]], max_deg)
  else
    return zero(u)
  end
end

# TODO: There could be cancellations un the weights of the normal bundle, e.g.
# 2e1 * e2 // e1 ** 2e2
# that don't carry over to cancel on the level of sinch.
# Thus, we really need a sinch factor for each weight of the normal bundle.
# Reimplement it that way.
function gkm_5d_free_strip_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64, res)

  @req !iszero(beta) "beta must not be zero"

  max_deg = 2*max_genus - 2
  
  # only nonzero exponents should be one
  b2 = rank(parent(beta))
  any(i -> beta[i] > 1, 1:b2) && return zero(u)

  
  # nonzero exponents should be adjacent
  first_nonzero = findfirst(i -> !iszero(beta[i]), 1:b2)
  last_nonzero = findlast(i -> !iszero(beta[i]), 1:b2)
  any(i -> iszero(beta[i]), first_nonzero+1:last_nonzero-1) && return zero(u)

  # factor the genus zero part
  g0 = evaluate(u^2 * res, vcat(t, [zero(u)]))
  # println("g0=$g0")
  num_fac = factor(numerator(g0))
  denom_fac = factor(denominator(g0))

  nfs = Vector{}()
  dfs = Vector{}()

  for (nf, e) in num_fac
    append!(nfs, repeat([u*nf], e))
  end
  for (df, e) in denom_fac
    append!(dfs, repeat([u*df], e))
  end
  
  # println("nfs = $nfs")
  # println("dfs = $dfs")
  # println("max_deg = $max_deg")

  return one(u) * unit(num_fac) * unit(denom_fac) * t_sq_bracket_frac(nfs, dfs, max_deg)

end