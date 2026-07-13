# run 1 implemented.
using Oscar, GKMtools, CY5

G = gkm_5d_P3_22()
b = curve_class(G, Edge(1, 2))

# Has this finished?

get_Omega_beta(G, [4*b], 3; check_predictions=true)