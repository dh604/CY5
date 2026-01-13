# Test GW vanishing for (O(-1,-1)+O(-1,-1) on P1xP1) x C

P1_P1 = hirzebruch_surface(NormalToricVariety, 0)
P = picard_group(P1_P1)
L = toric_line_bundle(P1_P1, -P[1] - P[2]) # O(-1,-1)
CY4 = total_space(gkm_vector_bundle_of_toric([L, L]))
G = CY5_from_CY4(CY4; equiCY=true)

b = curve_class(G, Edge(1, 2))
c = curve_class(G, Edge(1, 4))

# perform tests

get_Omega_beta(G, [4*b], 2; check_predictions = true)
# all zero

get_Omega_beta(G, [3*b], 2; check_predictions = true)
# all zero

get_Omega_beta(G, [2*b+2*c], 2; check_predictions = true)
# all zero

get_Omega_beta(G, [2*b+c], 2; check_predictions = true)
# all zero

get_Omega_beta(G, [3*b+c], 2; check_predictions = true)
# all zero