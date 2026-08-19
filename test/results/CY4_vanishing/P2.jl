# Test GW vanishing for (O(-1)+O(-2) on P^2) x C^1

using Oscar, GKMtools, GW_CY5

CY4 = total_space(vector_bundle_O(2, [-1, -2]))
G = CY5_from_CY4(CY4; equiCY=true)

b = curve_class(G, Edge(1, 2))

# perform tests

get_Omega_beta(G, [4*b], 3; check_predictions = true)
# all zero

get_Omega_beta(G, [3*b], 3; check_predictions = true)
# all zero