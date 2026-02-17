
@doc raw"""
    gkm_5d_p1_chain(steps::Vector{Vector{Int64}}, equiCY::Bool=false) -> GKMtools.AbstractGKM_graph

Return a 5d chain of $\mathbb{P}^1$s, where the degrees of the linear summands of the normal bundle
 along each $\mathbb{P}^1$ are a permutation of
`[2, -2, 0, 0, 0]` or `[2, -1, -1, 0, 0]`.

## Arguments
* `steps`: the vector containing for each $\mathbb{P}^1$ in the chain the degrees of its normal bundle splitting.
    The edges in the resulting GKM graph correspond to the places where 2 appears in the above splitting.
* `equiCY`: If this optional argument with default value `false` is set to `true` then the result will be linearized
    to be equivariantly Calabi-Yau.

!!! note
    The first step must have `2` in the last position.
    This ensures that the standalone flags at vertex `1` have weights $\epsilon_1,\dots,\epsilon_4$.
"""
function gkm_5d_p1_chain(steps::Vector{Vector{Int64}}; equiCY::Bool=false)
  @req all(step -> sort(step) in [[-2, 0, 0, 0, 2], [-1, -1, 0, 0, 2]], steps) "Step must be permutation of [2, -2, 0, 0, 0] or [-1, -1, 2, 0, 0]"
  
  l = length(steps)
  @req l>0 "Need at least one step."

  twos = [findfirst(s -> s == 2, step) for step in steps]
  @req all(i -> twos[i] != twos[i+1], 1:l-1) "2 must be in different positions for consecutive steps."
  @req twos[1] == 5 "The first step must have 2 in the last position"

  # Create a free module of rank 5
  M = free_module(ZZ, 5)
  g = gens(M)

  # Compute weights at each vertex
  n_verts = l + 1
  weights = Vector{Vector{typeof(g[1])}}(undef, n_verts)
  weights[1] = [g[1], g[2], g[3], g[4], g[5]]

  for i in 1:l
    w = weights[i]
    edge_w = w[twos[i]]
    weights[i+1] = [w[j] - steps[i][j] * edge_w for j in 1:5]
  end

  # Create GKM graph with all flags as standalone
  G = flags_only_gkm_graph(["$i" for i in 1:n_verts], M, weights)

  # Connect flags to form edges
  for i in 1:l
    connect_flags!(G, i, i+1, twos[i], twos[i])
  end

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

  set_attribute!(G, :example_type, :P1_chain_5d)
  return G
end