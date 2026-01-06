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
  uad = Vector{Int64}()
  for x in ups_and_downs
    if x > 0
      append!(uad, ones(Int64, x))
    else
      append!(uad, -ones(Int64, -x))
    end
  end

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
  f = ModuleHomomorphism(M, M, [g[2] + g[3], g[1], g[2]])
  G = substitute_torus(G, f)

  if equiCY
    M2 = free_module(ZZ, 2)
    g = gens(M2)
    f = ModuleHomomorphism(M, M2, [-g[1], g[2], -g[2]])
    G = substitute_torus(G, f)
  end

  return G
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