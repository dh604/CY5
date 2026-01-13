
@doc raw"""
    get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)

Return the conjectural prediction of $\Omega_\beta$, using the attribute `:example_type` of `G`.
"""
function get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)
  @req has_attribute(G, :example_type) "G has no prediction for Omega."
  et = get_attribute(G, :example_type)

  if et == :gkm_5d_strip_from_3d_CY
    return gkm_5d_strip_prediction(G, t, u, beta, max_genus)
  elseif et == :X_times_Ar
    return zero(u)
  elseif et == :closed_vertex
    return gkm_5d_closed_vertex_prediction(G, t, u, beta, max_genus)
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