using Oscar, GKMtools, GW_CY5

G = gkm_5d_P3_04()
b = curve_class(G, Edge(1, 2))


get_Omega_beta(G, [b, 2*b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [b, 3*b], 2; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [b], 5; check_predictions=true)
# All predictions hold.