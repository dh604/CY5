function gkm_3d_closed_vertex()

  ray_gens = [[1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1], [1, 2, 0], [1, 0, 2]]
  max_cones_list = Vector{Int64}[[1, 2, 3], [2, 3, 4], [2, 4, 5], [3, 4, 6]]
  N = normal_toric_variety(incidence_matrix(max_cones_list), ray_gens)
  G = gkm_graph_of_toric(N; small_torus=true)
  G.labels = ["1", "0", "2", "3"]

  # need QQ weights for substitution.
  G = convert_weights(G)
  
  # Substitute the torus to agree with Yannik's notation
  M = G.M
  g = gens(M)
  f = ModuleHomomorphism(M, M, [QQ(1)//QQ(2) * (-g[1] - g[2] - g[3]), -g[2], -g[3]])
  G = substitute_torus(G, f)

  return G
end

function gkm_5d_closed_vertex(;equiCY::Bool=false)

  CV3 = gkm_3d_closed_vertex()

  C2 = empty_gkm_graph(1, 2, ["1"])
  g = gens(C2.M)
  add_standalone_flag!(C2, 1, g[1])
  add_standalone_flag!(C2, 1, g[2])
  res = *(CV3, C2; calculateCurveClasses=false)
  if equiCY
    M = res.M
    g = gens(M)
    f = ModuleHomomorphism(M, M, [g[1], g[2], -g[1] - g[2] - 2*g[4] - 2*g[5], g[4], g[5]])
    res = substitute_torus(res, f)
    t = gens(res.equivariantCohomology.coeffRing)
    set_attribute!(res, :equiCY_substitution, [t[1], t[2], -t[1] - t[2] -2*(t[4] + t[5]), t[4], t[5]])
  end

  set_attribute!(res, :example_type, :closed_vertex)

  return res
end