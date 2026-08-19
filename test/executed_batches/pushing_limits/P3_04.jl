# run 1 implemented.
using Oscar, GKMtools, GW_CY5

G = gkm_5d_P3_04()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 3*b], 2; check_predictions=true) # was genus 1 before.
# Up to genus 2, l predictions hold.


# Old up-to-genus-1 prediction:

# Prediction holds for (1)
# Prediction holds for (3)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (1) => -18
#   (3) => -183778
