using Oscar, GKMtools, CY5

G = gkm_5d_strip([2, -1, 1]; equiCY=true)
b1 = curve_class(G, Edge(1, 2))
b2 = curve_class(G, Edge(2, 3))
b3 = curve_class(G, Edge(3, 4))

max_genus = 3

get_Omega_beta(G, [4*b1], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [4*b2], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [4*b3], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [4*b1+2*b2+4*b3], max_genus; check_predictions=true)

get_Omega_beta(G, [4*b1+4*b2+2*b3], max_genus; check_predictions=true)

get_Omega_beta(G, [4*b1+4*b2+4*b3], max_genus; check_predictions=true)