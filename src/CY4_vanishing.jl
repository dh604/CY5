@doc raw"""
    CY5_from_CY4(CY4::AbstractGKM_graph; equiCY::Bool=false)

Return $X\times\mathbb{C}^1$.

## Arguments
* `X`: The GKM graph of a 4-fold.
* `equiCY`: If this optional argument with default value `false` is set to `true`, the returned space will be linearized to be equivariantly Calabi-Yau.

!!! warning
    For `equiCY=true`, this function requires `X` to be Calabi-Yau, although not necessarily equivariantly.

# Example
Let us consider the example with $X=\text{Tot}(\mathcal{O}_{\mathbb{P}^2}(-1)\oplus \mathcal{O}_{\mathbb{P}^2}(-2))$.
```jldoctest
julia> X = total_space(vector_bundle_O(2, [-1, -2]))
GKM graph with 3 nodes, valency 4 and axial function:
2 -> 1 => (-1, 1, 0, 0, 0)
3 -> 1 => (-1, 0, 1, 0, 0)
3 -> 2 => (0, -1, 1, 0, 0)
Standalone flags:
1.3 => (0, 0, 0, 1, 0)
1.4 => (0, 0, 0, 0, 1)
2.3 => (1, -1, 0, 1, 0)
2.4 => (2, -2, 0, 0, 1)
3.3 => (1, 0, -1, 1, 0)
3.4 => (2, 0, -2, 0, 1)

julia> CY5_from_CY4(X)
GKM graph with 3 nodes, valency 5 and axial function:
2,1 -> 1,1 => (-1, 1, 0, 0, 0, 0)
3,1 -> 1,1 => (-1, 0, 1, 0, 0, 0)
3,1 -> 2,1 => (0, -1, 1, 0, 0, 0)
Standalone flags:
1,1.3 => (0, 0, 0, 1, 0, 0)
1,1.4 => (0, 0, 0, 0, 1, 0)
1,1.5 => (0, 0, 0, 0, 0, 1)
2,1.3 => (1, -1, 0, 1, 0, 0)
2,1.4 => (2, -2, 0, 0, 1, 0)
2,1.5 => (0, 0, 0, 0, 0, 1)
3,1.3 => (1, 0, -1, 1, 0, 0)
3,1.4 => (2, 0, -2, 0, 1, 0)
3,1.5 => (0, 0, 0, 0, 0, 1)

julia> CY5_from_CY4(X; equiCY=true) # Let us also see the equivariantly CY substitution.
GKM graph with 3 nodes, valency 5 and axial function:
2,1 -> 1,1 => (-1, 1, 0, 0, 0, 0)
3,1 -> 1,1 => (-1, 0, 1, 0, 0, 0)
3,1 -> 2,1 => (0, -1, 1, 0, 0, 0)
Standalone flags:
1,1.3 => (0, 0, 0, 1, 0, 0)
1,1.4 => (0, 0, 0, 0, 1, 0)
1,1.5 => (-2, 1, 1, -1, -1, 0)
2,1.3 => (1, -1, 0, 1, 0, 0)
2,1.4 => (2, -2, 0, 0, 1, 0)
2,1.5 => (-2, 1, 1, -1, -1, 0)
3,1.3 => (1, 0, -1, 1, 0, 0)
3,1.4 => (2, 0, -2, 0, 1, 0)
3,1.5 => (-2, 1, 1, -1, -1, 0)
```
"""
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