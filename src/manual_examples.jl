#
# TODO: fill in correct linearizations in this file, export functions, and document.
#

# TODO: document, noting that it is always CY5.
function gkm_5d_P2_111()

  G = total_space(vector_bundle_O(2, [-1, -1, -1]))
  M = G.M
  M2 = free_module(ZZ, 4)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [zero(M2), 2*g[1], 2*g[2], 2*g[3], 2*g[4], -2 * (g[1] + g[2] + g[3] + g[4])]) # TODO: find right linearization.
  G = substitute_torus(G, f)
  
  set_attribute!(res, :example_type, :P2_111)
  return G
end

# TODO: document, noting that it is always CY5.
function gkm_5d_P3_04()

  G = total_space(vector_bundle_O(3, [0, -4]))
  M = G.M
  M2 = free_module(ZZ, 4)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [zero(M2), 2*g[1], 2*g[2], 2*g[3], 2*g[4], -2 * (g[1] + g[2] + g[3] + g[4])]) # TODO: find right linearization.
  G = substitute_torus(G, f)
  
  set_attribute!(res, :example_type, :P3_04)
  return G
end

# TODO: document, noting that it is always CY5.
function gkm_5d_P3_13()

  G = total_space(vector_bundle_O(3, [-1, -3]))
  M = G.M
  M2 = free_module(ZZ, 4)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [zero(M2), 2*g[1], 2*g[2], 2*g[3], 2*g[4], -2 * (g[1] + g[2] + g[3] + g[4])]) # TODO: find right linearization.
  G = substitute_torus(G, f)
  
  set_attribute!(res, :example_type, :P3_13)
  return G
end

# TODO: document, noting that it is always CY5.
function gkm_5d_P3_22()

  G = total_space(vector_bundle_O(3, [-2, -2]))
  M = G.M
  M2 = free_module(ZZ, 4)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [zero(M2), 2*g[1], 2*g[2], 2*g[3], 2*g[4], -2 * (g[1] + g[2] + g[3] + g[4])]) # TODO: find right linearization.
  G = substitute_torus(G, f)
  
  set_attribute!(res, :example_type, :P3_22)
  return G
end