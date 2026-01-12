
@doc raw"""
    get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)

Return the conjectural prediction of $\Omega_\beta$, using the attribute `:example_type` of `G`.
"""
function get_Omega_prediction(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)
  @req has_attribute(G, :example_type) "G has no prediction for Omega."
  et = get_attribute(G, :example_type)

  if et == :gkm_5d_strip_from_3d_CY
    return gkm_5d_strip_prediction(G, t, u, beta, max_genus)
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


#TODO understand the sign difference in the following example:

# julia> G = gkm_5d_strip([-1, 1]; equiCY=true)
# GKM graph with 2 nodes, valency 5 and axial function:
# 2,1 -> 1,1 => (1, -1, 0, 0, 0)
# Standalone flags:
# 1,1.2 => (0, -1, 0, -1, -1)
# 1,1.3 => (1, 0, 0, 0, 0)
# 1,1.4 => (0, 0, 0, 1, 0)
# 1,1.5 => (0, 0, 0, 0, 1)
# 2,1.2 => (-1, 0, 0, -1, -1)
# 2,1.3 => (0, 1, 0, 0, 0)
# 2,1.4 => (0, 0, 0, 1, 0)
# 2,1.5 => (0, 0, 0, 0, 1)

# julia> b = curve_class(G, Edge(1, 2))
# (1)

# julia> get_Omega_beta(G, [b], 2)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => (7//5760*t4^4*u^4 + 1//576*t4^2*t5^2*u^4 - 1//24*t4^2*u^2 + 7//5760*t5^4*u^4 - 1//24*t5^2*u^2 + 1)//(t4*t5*u^2)

# julia> CY5.gkm_5d_strip_prediction_I_0(t, u, 3)
# (-7//5760*t4^4*u^4 + 1//576*t4^2*t5^2*u^4 - 1//24*t4^2*u^2 - 7//5760*t5^4*u^4 - 1//24*t5^2*u^2 + 1)//(t4*t5*u^2)