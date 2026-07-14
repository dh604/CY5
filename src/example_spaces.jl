# Example space from the article.

@doc raw"""
    gkm_graph_of_example_1_1()

Return the GKM graph of $\text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times\mathbb{C}^3$ as defined in [main_paper; Example 1.1](@cite).

# Example
```jldoctest
julia> G = gkm_graph_of_example_1_1()
GKM graph with 2 nodes, valency 5 and axial function:
2 -> 1 => (1, -1, 0, 0, 0, 0)
Standalone flags:
1.2 => (2, 0, 1, 0, 0, 0)
1.3 => (0, 0, 0, 1, 0, 0)
1.4 => (0, 0, 0, 0, 1, 0)
1.5 => (0, 0, 0, 0, 0, 1)
2.2 => (0, 2, 1, 0, 0, 0)
2.3 => (0, 0, 0, 1, 0, 0)
2.4 => (0, 0, 0, 0, 1, 0)
2.5 => (0, 0, 0, 0, 0, 1)
```
"""
function gkm_graph_of_example_1_1()
  M = free_module(ZZ, 6)
  g0, g1, g2, g3, g4, g5 = gens(M)
  vertex_labels = ["1", "2"]
  weights_at_0 = [g1-g0, g2+2*g0, g3, g4, g5]
  weights_at_1 = [g0-g1, g2+2*g1, g3, g4, g5]
  G = flags_only_gkm_graph(vertex_labels, M, [weights_at_0, weights_at_1])
  connect_flags!(G, 1, 2, 1, 1)
  return G
end