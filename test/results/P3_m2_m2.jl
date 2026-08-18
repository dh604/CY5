# Test independence of antidiagonal scaling parameter for O(-2)+O(-2) on P3.

using Oscar, GKMtools, GW_CY5

G = total_space(vector_bundle_O(3, [-2, -2]))
M = G.M
M2 = free_module(ZZ, 4)
g = gens(M2)
f = ModuleHomomorphism(M, M2, [zero(M2), -g[1], -g[2], -g[3], g[4], -g[4]])
G = substitute_torus(G, f)

# curve classes of A_r
b = curve_class(G, Edge(1, 2))

function is_indep_of_t4(x)::Bool
  x = x // 1
  !all(t -> iszero(div(t, gens(parent(t))[4])), terms(numerator(x))) && return false
  !all(t -> iszero(div(t, gens(parent(t))[4])), terms(denominator(x))) && return false
  return true
end

# run tests
D = get_Omega_beta(G, [b, 4*b, 3*b], 2)

is_indep_of_t4(D[b])
# true

is_indep_of_t4(D[2*b])
# true

is_indep_of_t4(D[3*b])
# false

# julia> D[b]
# 1//36*t1^2*u^2 + 1//18*t1*t2*u^2 + 1//18*t1*t3*u^2 + 1//36*t2^2*u^2 + 1//18*t2*t3*u^2 + 1//36*t3^2*u^2

# julia> D[2*b]
# -107//90*t1^2*u^2 - 107//45*t1*t2*u^2 - 107//45*t1*t3*u^2 - 107//90*t2^2*u^2 - 107//45*t2*t3*u^2 - 107//90*t3^2*u^2 + 2

# julia> D[3*b]
# 1919//12*t1^2*u^2 + 1439//6*t1*t2*u^2 + 1439//6*t1*t3*u^2 + 1919//12*t2^2*u^2 + 1439//6*t2*t3*u^2 + 1919//12*t3^2*u^2 - 40*t4^2*u^2 + 32