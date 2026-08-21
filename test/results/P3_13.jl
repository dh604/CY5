using Oscar, GKMtools, GW_CY5

G = gkm_5d_P3_13()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b], 5; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*b], 5; check_predictions=true)
# All predictions hold.