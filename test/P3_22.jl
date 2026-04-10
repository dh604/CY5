G = gkm_5d_P3_22()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 2*b, 3*b], 1; check_predictions=true)

# Prediction holds for (2)
# Prediction holds for (1)
# Prediction holds for (3)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 3 entries:
#   (2) => 2
#   (1) => 0
#   (3) => 32

get_Omega_beta(G, [b, 2*b], 3; check_predictions=true)

# Prediction holds for (2)
# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (2) => 2
#   (1) => 0

get_Omega_beta(G, [b], 5; check_predictions=true)

# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => 0