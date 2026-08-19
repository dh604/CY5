# Old first tests for O(-1)+O(-1)+O(-1) on P2.

using Oscar, GKMtools, GW_CY5

G = total_space(vector_bundle_O(2, [-1, -1, -1]))
M = G.M
M2 = free_module(ZZ, 4)
g = gens(M2)
f = ModuleHomomorphism(M, M2, [zero(M2), -g[1], -g[2], g[3], g[4], -g[1] - g[2] - g[3] - g[4]])
G = substitute_torus(G, f)

b = curve_class(G, Edge(1, 2))

# function is_indep_of_t4(x)::Bool
#   x = x // 1
#   !all(t -> iszero(div(t, gens(parent(t))[4])), terms(numerator(x))) && return false
#   !all(t -> iszero(div(t, gens(parent(t))[4])), terms(denominator(x))) && return false
#   return true
# end

# run tests
D = get_Omega_beta(G, [b, 4*b, 3*b ,5*b], 2)

# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 4 entries:
#   (1) => 1//128*t1^2*u^2 - 1//128*t1*t2*u^2 + 1//128*t2^2*u^2 - 1//8
#   (2) => -1//32*t1^2*u^2 + 7//128*t1*t2*u^2 + 3//128*t1*t3*u^2 + 3//128*t1*t4*u^2 - 1//32*t2^2*u^2 + 3//128*t2*t3*u^2 + 3//128*t2*t4*u^2 + 3//128*t3^2*u^2 + 3//128*t3*t4*u^2 + 3//128*t4^2*u^2 + 1//8
#   (3) => -9//128*t1*t2*u^2 - 9//128*t1*t3*u^2 - 9//128*t1*t4*u^2 - 9//128*t2*t3*u^2 - 9//128*t2*t4*u^2 - 9//128*t3^2*u^2 - 9//128*t3*t4*u^2 - 9//128*t4^2*u^2
#   (4) => 3//32*t1*t2*u^2 + 3//32*t1*t3*u^2 + 3//32*t1*t4*u^2 + 3//32*t2*t3*u^2 + 3//32*t2*t4*u^2 + 3//32*t3^2*u^2 + 3//32*t3*t4*u^2 + 3//32*t4^2*u^2
#   (5) => -15//64*t1*t2*u^2 - 15//64*t1*t3*u^2 - 15//64*t1*t4*u^2 - 15//64*t2*t3*u^2 - 15//64*t2*t4*u^2 - 15//64*t3^2*u^2 - 15//64*t3*t4*u^2 - 15//64*t4^2*u^2
#   (6) => 27//64*t1*t2*u^2 + 27//64*t1*t3*u^2 + 27//64*t1*t4*u^2 + 27//64*t2*t3*u^2 + 27//64*t2*t4*u^2 + 27//64*t3^2*u^2 + 27//64*t3*t4*u^2 + 27//64*t4^2*u^2