@doc raw"""
    gkm_3d_strip(ups_and_downs::Vector{Int64}; equiCY::Bool=false)::GKMtools.AbstractGKM_graph

Create the GKM graph of the strip geometry described by `ups_and_downs`.
This means that there will be `ups_and_downs[1]` flags up (if positive) or down (if negative),
followed by `ups_and_downs[2]` flags up or down, and so on.

# Weight convention:
Flags pointing up have weight $\epsilon_2$, flags pointing down have weight $\epsilon_3$, and
the flag at vertex 1 pointing left has weight $\epsilon_1$.

If the optional argument `equiCY` is set to `true`, then the weights are as in (4) of the draft.
"""
function gkm_3d_strip(ups_and_downs::Vector{Int64}; equiCY::Bool=false)::GKMtools.AbstractGKM_graph

  @req all(i -> !iszero(i), ups_and_downs) "ups_and_downs must not contain a zero entry"
  @req length(ups_and_downs) >= 1 "ups_and_downs must not be empty"

  # convert into vector of 1 and -1
  uad = uad_only_ones(ups_and_downs)

  ray_gens = [[1, 0, 0], [1, 0, 1]]
  max_cones_list = Vector{Vector{Int64}}()
  last_up_index = 1
  last_down_index = 2
  last_up_coord = 0
  last_down_coord = 0
  for x in uad
    if x == 1 # up
      push!(ray_gens, [1, last_up_coord + 1, 0])
      push!(max_cones_list, [last_up_index, last_down_index, length(ray_gens)])
      last_up_coord += 1
      last_up_index = length(ray_gens)
    elseif x == -1 # down
      push!(ray_gens, [1, last_down_coord + 1, 1])
      push!(max_cones_list, [last_up_index, last_down_index, length(ray_gens)])
      last_down_coord += 1
      last_down_index = length(ray_gens)
    end
  end
  # println("ray_gens: $ray_gens")
  # println("max_cones_list: $max_cones_list")
  N = normal_toric_variety(incidence_matrix(max_cones_list), ray_gens)
  G = gkm_graph_of_toric(N; small_torus=true)
  
  # Substitute the torus to agree with Yannik's notation
  M = G.M
  g = gens(M)
  # f = ModuleHomomorphism(M, M, [g[2] + g[3], g[1], g[2]])
  f = ModuleHomomorphism(M, M, [-g[2] - g[3], -g[1], -g[2]])
  G = substitute_torus(G, f)

  if equiCY
    M2 = free_module(ZZ, 2)
    g = gens(M2)
    f = ModuleHomomorphism(M, M2, [-g[1], g[2], -g[2]])
    G = substitute_torus(G, f)
  end

  return G
end


function uad_only_ones(ups_and_downs::Vector{Int64})::Vector{Int64}
  uad = Vector{Int64}()
  for x in ups_and_downs
    if x > 0
      append!(uad, ones(Int64, x))
    else
      append!(uad, -ones(Int64, -x))
    end
  end
  return uad
end


function test_strips(l::Int64)
  r = abs.(rand(Int64, l)) .% 9 .+ 1
  r2 = abs.(rand(Int64, l)) .% 2
  r .*= (-1).^r2
  G = gkm_3d_strip(r)
  @req isvalid(G) "G is invalid for r=$r"
  @req is_gkm_class(chern_class(G, 1), G) "c1 not GKM for r=$r"
  print_curve_classes(G)
end

@doc raw"""
    gkm_5d_strip(ups_and_downs::Vector{Int64}; equiCY::Bool=false)::GKMtools.AbstractGKM_graph

Return the product of the output of `gkm_3d_strip` with $\mathbb{C}^2$, where the trivial factor
has weights $\epsilon_4, \epsilon_5$.
"""
function gkm_5d_strip(ups_and_downs::Vector{Int64}; equiCY::Bool=false)::GKMtools.AbstractGKM_graph
  three_d_strip = gkm_3d_strip(ups_and_downs)
  C2 = empty_gkm_graph(1, 2, ["1"])
  g = gens(C2.M)
  add_standalone_flag!(C2, 1, g[1])
  add_standalone_flag!(C2, 1, g[2])
  res = *(three_d_strip, C2; calculateCurveClasses=false)
  if equiCY
    M = res.M
    g = gens(M)
    # f = ModuleHomomorphism(M, M, [g[1], g[2], g[3], g[4], -g[2]-g[3]-g[4]])
    f = ModuleHomomorphism(M, M, [g[1], g[2], -g[2]-g[4]-g[5], g[4], g[5]])
    res = substitute_torus(res, f)
    t = gens(res.equivariantCohomology.coeffRing)
    set_attribute!(res, :equiCY_substitution, [t[1], t[2], -t[2]-t[4]-t[5], t[4], t[5]])
  end

  # Save I_u, I_d, and I_0 information
  # if this throws an error, it means that the toric part of the code didn't label the vertices 1, 2, 3, ... along the strip.
  strip_dict = Dict{CC, Symbol}()
  nv = n_vertices(res.g)
  uad = uad_only_ones(ups_and_downs)
  for i in 1:nv
    for j in i+1:nv
      b = sum(k -> curve_class(res, Edge(k, k+1)), i:j-1)
      if uad[i] != uad[j]
        strip_dict[b] = :I_0
      elseif uad[i] == 1
        strip_dict[b] = :I_u
      else # uad[i] = uad[j] = -1
        strip_dict[b] = :I_d
      end
    end
  end
  set_attribute!(res, :prediction_data, strip_dict)
  set_attribute!(res, :example_type, :gkm_5d_strip_from_3d_CY)
  return res
end

# # TODO: implement
# function gkm_p1_chain()

# end