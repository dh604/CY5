
@doc raw"""
    gkm_3d_gauge(N::Int64, m::Int64)

Return the GKM graph of the CY3 space $X_{N,m}$ from [main_paper; Section 5](@cite).
The returned space is Calabi-Yau but not equivariantly Calabi-Yau.

This function is used by [`gkm_5d_gauge`](@ref) to construct the GKM graph of $\mathcal{Z}_{N,m} = X_{N,m}\times\mathbb{C}^2$.

# Example
```jldoctest
julia> G = gkm_3d_gauge(3, 5)
GKM graph with 6 nodes, valency 3 and axial function:
v2 -> v1 => (1, 1, 0)
v3 -> v2 => (1, 2, 0)
w1 -> v1 => (0, 0, -1)
w2 -> v2 => (0, 0, -1)
w2 -> w1 => (1, 1, 6)
w3 -> v3 => (0, 0, -1)
w3 -> w2 => (1, 2, 4)
Standalone flags:
v1.2 => (1, 0, 0)
v3.3 => (-1, -3, 0)
w1.2 => (1, 0, 8)
w3.3 => (-1, -3, -2)

julia> print_curve_classes(G)
v2 -> v1: (1, 0, 0), Chern number: 0
v3 -> v2: (0, 0, 1), Chern number: 0
w1 -> v1: (6, 1, 4), Chern number: 0
w2 -> v2: (0, 1, 4), Chern number: 0
w2 -> w1: (1, 0, 0), Chern number: 0
w3 -> v3: (0, 1, 0), Chern number: 0
w3 -> w2: (0, 0, 1), Chern number: 0

julia> first_chern_class(G)
(-t2 + t3)*e[1] + (-t2 + t3)*e[2] + (-t2 + t3)*e[3] + (-t2 + t3)*e[4] + (-t2 + t3)*e[5] + (-t2 + t3)*e[6]
```
"""
function gkm_3d_gauge(N::Int64, m::Int64)

  @req N > 0 "N must be positive"

  labels = vcat(["v$i" for i in 1:N], ["w$i" for i in 1:N])
  G = empty_gkm_graph(2*N, 3, labels)
  g1, g2, g3 = gens(G.M)

  # Add edges labeled Q_i
  for i in 1:N
    add_edge!(G, "v$i", "w$i", g3)
  end

  # Add standalone flags on top
  add_standalone_flag!(G, "v1", g1)
  add_standalone_flag!(G, "w1", g1 + (N+m)*g3)

  # Add edges labeled \tilde{Q}_i
  for i in 1:N-1
    add_edge!(G, "v$(i+1)", "v$i", g1 + i*g2)
    add_edge!(G, "w$(i+1)", "w$i", g1 + i*g2 + (N+m-2*i)*g3)
  end

  # Add standalone flags on bottom
  add_standalone_flag!(G, "v$N", -(g1 + N*g2))
  add_standalone_flag!(G, "w$N", -(g1 + N*g2 + (-N+m)*g3))

  return G
end


@doc raw"""
    gkm_5d_gauge(N::Int64, m::Int64; equiCY::Bool=false)

Return the GKM graph of the CY5 space $\mathcal{Z}_{N,m}$ from [main_paper; Section 5](@cite).
If the optional argument `equiCY` is set to `true` then the substitution $\epsilon_2=\epsilon_3+\epsilon_4+\epsilon_5$ is applied to make
the result equivariantly Calabi-Yau.

# Examples
Let us first see the non-equivariantly Calabi-Yau example.
```jldoctest
julia> G = gkm_5d_gauge(3, 5)
GKM graph with 6 nodes, valency 5 and axial function:
v2 -> v1 => (1, 1, 0, 0, 0)
v3 -> v2 => (1, 2, 0, 0, 0)
w1 -> v1 => (0, 0, -1, 0, 0)
w2 -> v2 => (0, 0, -1, 0, 0)
w2 -> w1 => (1, 1, 6, 0, 0)
w3 -> v3 => (0, 0, -1, 0, 0)
w3 -> w2 => (1, 2, 4, 0, 0)
Standalone flags:
v1.2 => (1, 0, 0, 0, 0)
v1.4 => (0, 0, 0, 1, 0)
v1.5 => (0, 0, 0, 0, 1)
v2.4 => (0, 0, 0, 1, 0)
v2.5 => (0, 0, 0, 0, 1)
v3.3 => (-1, -3, 0, 0, 0)
v3.4 => (0, 0, 0, 1, 0)
v3.5 => (0, 0, 0, 0, 1)
w1.2 => (1, 0, 8, 0, 0)
w1.4 => (0, 0, 0, 1, 0)
w1.5 => (0, 0, 0, 0, 1)
w2.4 => (0, 0, 0, 1, 0)
w2.5 => (0, 0, 0, 0, 1)
w3.3 => (-1, -3, -2, 0, 0)
w3.4 => (0, 0, 0, 1, 0)
w3.5 => (0, 0, 0, 0, 1)

julia> first_chern_class(G)
(-t2 + t3 + t4 + t5)*e[1] + (-t2 + t3 + t4 + t5)*e[2] + (-t2 + t3 + t4 + t5)*e[3] + (-t2 + t3 + t4 + t5)*e[4] + (-t2 + t3 + t4 + t5)*e[5] + (-t2 + t3 + t4 + t5)*e[6]
```
Second, let us see the equivariantly Calabi-Yau example.
```jldoctest
julia> G = gkm_5d_gauge(3, 5; equiCY=true)
GKM graph with 6 nodes, valency 5 and axial function:
v2 -> v1 => (1, 0, 1, 1, 1)
v3 -> v2 => (1, 0, 2, 2, 2)
w1 -> v1 => (0, 0, -1, 0, 0)
w2 -> v2 => (0, 0, -1, 0, 0)
w2 -> w1 => (1, 0, 7, 1, 1)
w3 -> v3 => (0, 0, -1, 0, 0)
w3 -> w2 => (1, 0, 6, 2, 2)
Standalone flags:
v1.2 => (1, 0, 0, 0, 0)
v1.4 => (0, 0, 0, 1, 0)
v1.5 => (0, 0, 0, 0, 1)
v2.4 => (0, 0, 0, 1, 0)
v2.5 => (0, 0, 0, 0, 1)
v3.3 => (-1, 0, -3, -3, -3)
v3.4 => (0, 0, 0, 1, 0)
v3.5 => (0, 0, 0, 0, 1)
w1.2 => (1, 0, 8, 0, 0)
w1.4 => (0, 0, 0, 1, 0)
w1.5 => (0, 0, 0, 0, 1)
w2.4 => (0, 0, 0, 1, 0)
w2.5 => (0, 0, 0, 0, 1)
w3.3 => (-1, 0, -5, -3, -3)
w3.4 => (0, 0, 0, 1, 0)
w3.5 => (0, 0, 0, 0, 1)

julia> first_chern_class(G)
0
```
"""
function gkm_5d_gauge(N::Int64, m::Int64; equiCY::Bool=false)

  G3 = gkm_3d_gauge(N, m)

  C2 = empty_gkm_graph(1, 2, ["1"])
  g = gens(C2.M)
  add_standalone_flag!(C2, 1, g[1])
  add_standalone_flag!(C2, 1, g[2])
  res = *(G3, C2; calculateCurveClasses=false)
  res.labels = vcat(["v$i" for i in 1:N], ["w$i" for i in 1:N])
  if equiCY
    M = res.M
    g = gens(M)
    # Use substitution -e2 + e3 + e4 + e5 = 0 for e2.
    f = ModuleHomomorphism(M, M, [g[1], g[3]+g[4]+g[5], g[3], g[4], g[5]])
    res = substitute_torus(res, f)
    t = gens(res.equivariantCohomology.coeffRing)
    set_attribute!(res, :equiCY_substitution, [t[1], t[3]+t[4]+t[5], t[3], t[4], t[5]])
  end

  set_attribute!(res, :example_type, :gauge)
  set_attribute!(res, :N, N)
  set_attribute!(res, :m, m)

  return res
end