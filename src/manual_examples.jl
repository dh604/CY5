@doc raw"""
    gkm_5d_P2_111()

Return the GKM graph of $\text{Tot}\mathcal{O}_{\mathbb{P}^2}(-1)^{\oplus 3}$ with the specific equivariantly CY linearization as in the example below.

# Example

```jldoctest
julia> G = gkm_5d_P2_111()
GKM graph with 3 nodes, valency 5 and axial function:
2 -> 1 => (2, -2, 0, 0, 0)
3 -> 1 => (2, 0, -2, 0, 0)
3 -> 2 => (0, 2, -2, 0, 0)
Standalone flags:
1.3 => (2, 0, 0, 2, 0)
1.4 => (2, 0, 0, 0, 2)
1.5 => (0, -2, -2, -2, -2)
2.3 => (0, 2, 0, 2, 0)
2.4 => (0, 2, 0, 0, 2)
2.5 => (-2, 0, -2, -2, -2)
3.3 => (0, 0, 2, 2, 0)
3.4 => (0, 0, 2, 0, 2)
3.5 => (-2, -2, 0, -2, -2)

julia> first_chern_class(G)
0
```
"""
function gkm_5d_P2_111()

  G = total_space(vector_bundle_O(2, [-1, -1, -1]))
  M = G.M
  M2 = free_module(ZZ, 5)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [-2*g[1], -2*g[2], -2*g[3], 2*(g[4]+g[1]), 2*(g[5]+g[1]), -2*(g[2]+g[3]+g[4]+g[5])])
  G = substitute_torus(G, f)
  
  set_attribute!(G, :example_type, :P2_111)
  return G
end

@doc raw"""
    gkm_5d_P3_04()

Return the GKM graph of $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(0)\oplus \mathcal{O}_{\mathbb{P}^3}(-4))$ with the specific equivariantly CY linearization as in the example below.

# Example

```jldoctest
julia> G = gkm_5d_P3_04()
GKM graph with 4 nodes, valency 5 and axial function:
2 -> 1 => (2, -2, 0, 0, 0)
3 -> 1 => (2, 0, -2, 0, 0)
3 -> 2 => (0, 2, -2, 0, 0)
4 -> 1 => (2, 0, 0, -2, 0)
4 -> 2 => (0, 2, 0, -2, 0)
4 -> 3 => (0, 0, 2, -2, 0)
Standalone flags:
1.4 => (0, 0, 0, 0, 2)
1.5 => (6, -2, -2, -2, -2)
2.4 => (0, 0, 0, 0, 2)
2.5 => (-2, 6, -2, -2, -2)
3.4 => (0, 0, 0, 0, 2)
3.5 => (-2, -2, 6, -2, -2)
4.4 => (0, 0, 0, 0, 2)
4.5 => (-2, -2, -2, 6, -2)

julia> first_chern_class(G)
0
```
"""
function gkm_5d_P3_04()

  G = total_space(vector_bundle_O(3, [0, -4]))
  M = G.M
  M2 = free_module(ZZ, 5)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [-2*g[1], -2*g[2], -2*g[3], -2*g[4], 2*(g[5]), -2*(-3*g[1]+g[2]+g[3]+g[4]+g[5])])
  G = substitute_torus(G, f)
  
  set_attribute!(G, :example_type, :P3_04)
  return G
end

@doc raw"""
    gkm_5d_P3_13()

Return the GKM graph of $\text{Tot}(\mathcal{O}_{\mathbb{P}^3}(-1)\oplus \mathcal{O}_{\mathbb{P}^3}(-3))$ with the specific equivariantly CY linearization as in the example below.

# Example

```jldoctest
julia> G = gkm_5d_P3_13()
GKM graph with 4 nodes, valency 5 and axial function:
2 -> 1 => (2, -2, 0, 0, 0)
3 -> 1 => (2, 0, -2, 0, 0)
3 -> 2 => (0, 2, -2, 0, 0)
4 -> 1 => (2, 0, 0, -2, 0)
4 -> 2 => (0, 2, 0, -2, 0)
4 -> 3 => (0, 0, 2, -2, 0)
Standalone flags:
1.4 => (2, 0, 0, 0, 2)
1.5 => (4, -2, -2, -2, -2)
2.4 => (0, 2, 0, 0, 2)
2.5 => (-2, 4, -2, -2, -2)
3.4 => (0, 0, 2, 0, 2)
3.5 => (-2, -2, 4, -2, -2)
4.4 => (0, 0, 0, 2, 2)
4.5 => (-2, -2, -2, 4, -2)

julia> first_chern_class(G)
0
```
"""
function gkm_5d_P3_13()

  G = total_space(vector_bundle_O(3, [-1, -3]))
  M = G.M
  M2 = free_module(ZZ, 5)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [-2*g[1], -2*g[2], -2*g[3], -2*g[4], 2*(g[5]+g[1]), -2*(-2*g[1]+g[2]+g[3]+g[4]+g[5])])
  G = substitute_torus(G, f)
  
  set_attribute!(G, :example_type, :P3_13)
  return G
end

@doc raw"""
    gkm_5d_P3_22()

Return the GKM graph of $\text{Tot}\mathcal{O}_{\mathbb{P}^3}(-2)^{\oplus 2}$ with the specific equivariantly CY linearization as in the example below.

# Example

```jldoctest
julia> G = gkm_5d_P3_22()
GKM graph with 4 nodes, valency 5 and axial function:
2 -> 1 => (2, -2, 0, 0, 0)
3 -> 1 => (2, 0, -2, 0, 0)
3 -> 2 => (0, 2, -2, 0, 0)
4 -> 1 => (2, 0, 0, -2, 0)
4 -> 2 => (0, 2, 0, -2, 0)
4 -> 3 => (0, 0, 2, -2, 0)
Standalone flags:
1.4 => (4, 0, 0, 0, 2)
1.5 => (2, -2, -2, -2, -2)
2.4 => (0, 4, 0, 0, 2)
2.5 => (-2, 2, -2, -2, -2)
3.4 => (0, 0, 4, 0, 2)
3.5 => (-2, -2, 2, -2, -2)
4.4 => (0, 0, 0, 4, 2)
4.5 => (-2, -2, -2, 2, -2)

julia> first_chern_class(G)
0
```
"""
function gkm_5d_P3_22()

  G = total_space(vector_bundle_O(3, [-2, -2]))
  M = G.M
  M2 = free_module(ZZ, 5)
  g = gens(M2)
  f = ModuleHomomorphism(M, M2, [-2*g[1], -2*g[2], -2*g[3], -2*g[4], 2*(g[5]+2*g[1]), -2*(-1*g[1]+g[2]+g[3]+g[4]+g[5])])
  G = substitute_torus(G, f)
  
  set_attribute!(G, :example_type, :P3_22)
  return G
end