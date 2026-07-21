# run 1 implemented.
using Oscar, GKMtools, CY5

# Test conjecture for Tot(K_{Hirzebruch_5}) x A_r

r = 3

Y = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5))

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of H5
H1 = curve_class(G, "1,1", "2,1")
H2 = curve_class(G, "1,1", "4,1")

# run tests

get_Omega_beta(G, [3*b], 3; check_predictions=true)
# fails for b, holds for 3*b

get_Omega_beta(G, [H2 + b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*H2 + b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [3*H2 + b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [H2 + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*H2 + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [3*H2 + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [H2 + b + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*H2 + b + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [3*H2 + b + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [H2 + 2*b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*H2 + 2*b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [H2 + 2*c], 3; check_predictions=true)
# All predictions hold.