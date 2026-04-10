G = gkm_5d_P3_13()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 2*b], 3; check_predictions=true)

# Prediction holds for (2)
# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (2) => 13//4*t1^4*u^4 + 989//192*t1^3*t2*u^4 + 989//192*t1^3*t3*u^4 + 989//192*t1^3*t4*u^4 + 1821//64*t1^3*t5*u^4 + 413//64*t1^2*t2^2*u^4 - 345//64*t1^2*t2*t3*u^4 - 345//64*t1^2*t2*t4*u^4 + 1125//64*t1^2*t2*t5*u^4 + 413//64*t1^2*t3^…
#   (1) => -1//48*t1^3*t2*u^4 - 1//48*t1^3*t3*u^4 - 1//48*t1^3*t4*u^4 - 1//16*t1^3*t5*u^4 + 1//64*t1^2*t2^2*u^4 + 1//32*t1^2*t2*t3*u^4 + 1//32*t1^2*t2*t4*u^4 + 1//32*t1^2*t2*t5*u^4 + 1//64*t1^2*t3^2*u^4 + 1//32*t1^2*t3*t4*u^4 + 1//32*t1…


get_Omega_beta(G, [b], 5; check_predictions=true)

# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => -79//161280*t1^7*t2*u^8 - 79//161280*t1^7*t3*u^8 - 79//161280*t1^7*t4*u^8 - 79//53760*t1^7*t5*u^8 + 47//46080*t1^6*t2^2*u^8 + 47//23040*t1^6*t2*t3*u^8 + 47//23040*t1^6*t2*t4*u^8 + 31//11520*t1^6*t2*t5*u^8 + 47//46080*t1^6*t3^…
