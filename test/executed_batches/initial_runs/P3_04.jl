G = gkm_5d_P3_04()
b = curve_class(G, Edge(1, 2))


get_Omega_beta(G, [b, 2*b], 3; check_predictions=true)

# Prediction holds for (2)
# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (2) => -3576*t1^4*u^4 + 4768*t1^3*t2*u^4 + 4768*t1^3*t3*u^4 + 4768*t1^3*t4*u^4 - 11120*t1^2*t2^2*u^4 + 3968*t1^2*t2*t3*u^4 + 3968*t1^2*t2*t4*u^4 - 11120*t1^2*t3^2*u^4 + 3968*t1^2*t3*t4*u^4 - 11120*t1^2*t4^2*u^4 - 28560*t1^2*t5^2*u^4…
#   (1) => -3//4*t1^4*u^4 + t1^3*t2*u^4 + t1^3*t3*u^4 + t1^3*t4*u^4 - 9//2*t1^2*t2^2*u^4 + 3*t1^2*t2*t3*u^4 + 3*t1^2*t2*t4*u^4 - 9//2*t1^2*t3^2*u^4 + 3*t1^2*t3*t4*u^4 - 9//2*t1^2*t4^2*u^4 - 12*t1^2*t5^2*u^4 - 9*t1^2*u^2 + t1*t2^3*u^4 + …


get_Omega_beta(G, [b, 3*b], 1; check_predictions=true)

# Prediction holds for (1)
# Prediction holds for (3)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (1) => -18
#   (3) => -183778

get_Omega_beta(G, [b], 5; check_predictions=true)

# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => -1//2240*t1^8*u^8 + 1//840*t1^7*t2*u^8 + 1//840*t1^7*t3*u^8 + 1//840*t1^7*t4*u^8 - 1//80*t1^6*t2^2*u^8 + 1//120*t1^6*t2*t3*u^8 + 1//120*t1^6*t2*t4*u^8 - 1//80*t1^6*t3^2*u^8 + 1//120*t1^6*t3*t4*u^8 - 1//80*t1^6*t4^2*u^8 - 1//3…
