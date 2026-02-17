G = gkm_5d_p1_chain([[0, 0, -2, 0, 2], [0, 0, -2, 2, 0]]; equiCY=true)
b1 = curve_class(G, Edge(1, 2))
b2 = curve_class(G, Edge(2, 3))
# b3 = curve_class(G, Edge(3, 4))

get_Omega_beta(G, [b2, b1, 2*b2+2*b1], 2; check_predictions=true)