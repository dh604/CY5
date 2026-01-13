function CY5_from_CY4(CY4::AbstractGKM_graph; equiCY::Bool=false)

  C1 = empty_gkm_graph(1, 1, ["1"])
  g = gens(C1.M)
  add_standalone_flag!(C1, 1, g[1])
  G = *(CY4, C1; calculateCurveClasses=false)
  if equiCY
    M = G.M
    g = gens(M)
    c1 = sum(G.weights_at_vertex[1])
    subst_vect = [g[i] for i in 1:length(g)]
    subst_vect[length(g)] -= c1
    f = ModuleHomomorphism(M, M, subst_vect)
    G = substitute_torus(G, f)
    t = gens(G.equivariantCohomology.coeffRing)
    equiCY_subst = [t[i] for i in 1:length(g)]
    equiCY_subst[length(g)] -= sum(i -> c1[i] * t[i], 1:length(g))
    set_attribute!(G, :equiCY_substitution, equiCY_subst)
  end

  set_attribute!(G, :example_type, :CY5_from_CY4)

  return G
end