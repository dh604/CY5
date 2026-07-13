# run 1 implemented.
using Oscar, GKMtools, CY5

G = gkm_5d_P3_22()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 2*b, 3*b], 3; check_predictions=true) # beta=3b should be non-zero for large enough g (is g=3 large enough?) Yes. In fact, genus 1 is alreaday non-zero for 3*b.
# All predictions hold.

get_Omega_beta(G, [4*b], 2; check_predictions=true) # should be non-zero for some small g (is 2 large enough?)
# Crashed as we only have predictions up to 3*b.

# get_Omega_beta(G, [b, 2*b], 3; check_predictions=true)


# get_Omega_beta(G, [b], 5; check_predictions=true)

# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => 0