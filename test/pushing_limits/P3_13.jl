# run 1 implemented.
using Oscar, GKMtools, CY5

G = gkm_5d_P3_13()
b = curve_class(G, Edge(1, 2))

# ERRORED: only have prediction for b=1, 2.
get_Omega_beta(G, [b, 3*b], 2; check_predictions=true)