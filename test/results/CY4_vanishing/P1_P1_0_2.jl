# Test GW vanishing for (O(-2,0)+O(0,-2) on P1xP1) x C

using Oscar, GKMtools, CY5

P1_P1 = hirzebruch_surface(NormalToricVariety, 0)
P = picard_group(P1_P1)
L1 = toric_line_bundle(P1_P1, -2*P[1]) # O(-2, 0)
L2 = toric_line_bundle(P1_P1, -2*P[2]) # O(0, -2)
CY4 = total_space(gkm_vector_bundle_of_toric([L1, L2]))
G = CY5_from_CY4(CY4; equiCY=true)

b = curve_class(G, Edge(1, 2))
c = curve_class(G, Edge(1, 4))

# perform tests

get_Omega_beta(G, [4*b], 3; check_predictions = true)
# non-zero for b, zero for 2b, 4b

get_Omega_beta(G, [3*b], 3; check_predictions = true)
# nonzero for b, zero for 3b

get_Omega_beta(G, [2*b+2*c], 3; check_predictions = true)
# all zero

get_Omega_beta(G, [2*b+c], 3; check_predictions = true)
# all zero

get_Omega_beta(G, [3*b+c], 3; check_predictions = true)
# all zero
