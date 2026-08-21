# Julia code to test GW vanishing for O(-1)+O(-2)+O on P^2

using Oscar, GKMtools, GW_CY5

CY4 = total_space(vector_bundle_O(2, [-1, -2]))
G = CY5_from_CY4(CY4; equiCY=true)

b = curve_class(G, Edge(1, 2))
max_genus = 3

get_Omega_beta(G, [b, 2*b, 3*b, 4*b], max_genus; check_predictions = true)
# All predictions hold