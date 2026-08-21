using Oscar, GKMtools, GW_CY5

G = gkm_5d_P3_04()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 3*b], 3; check_predictions=true)
# All predictions hold.

# Prediction holds for (1)
# Prediction holds for (3)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (1) => -3//4*t1^4*u^4 + t1^3*t2*u^4 + t1^3*t3*u^4 + t1^3*t4*u^4 - 9//2*t1^2*t2^2*u^4 + 3*t1^2*t2*t3*u^4 + 3*t1^2*t2*t4*u^4 - 9//2*…
#   (3) => -12352203//4*t1^4*u^4 + 4117401*t1^3*t2*u^4 + 4117401*t1^3*t3*u^4 + 4117401*t1^3*t4*u^4 + 243936*t1^3*t5*u^4 - 17333601//2*…