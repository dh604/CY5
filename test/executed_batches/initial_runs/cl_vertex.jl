# Test the conjecture for the 5d closed vertex.

G = gkm_5d_closed_vertex(;equiCY=true)

b1 = curve_class(G, "0,1", "1,1")
b2 = curve_class(G, "0,1", "2,1")
b3 = curve_class(G, "0,1", "3,1")

get_Omega_beta(G, [4*b1, 6*b1], 2; check_predictions=true)
# all predictions hold

get_Omega_beta(G, [2*(b1+b2)], 2; check_predictions=true)
# all predictions hold

get_Omega_beta(G, [3*(b1+b2)], 2; check_predictions=true)
# all predictions hold

get_Omega_beta(G, [2*(b1+b2+b3)], 2; check_predictions=true)
# all predictions hold.

get_Omega_beta(G, [4*b1+2*b2], 2; check_predictions=true)
# all predictions hold.

get_Omega_beta(G, [b1+2*b2 + b3], 2; check_predictions=true)
# all predictions hold.

get_Omega_beta(G, [b1+3*b2 + b3], 2; check_predictions=true)
# all predictions hold.