using Oscar, GKMtools, CY5

# Test conjecture for Tot(K_{P2}) x A_r

r = 3

Y = projective_space(GKM_graph, 2)

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of P2
H = curve_class(G, "1,1", "2,1")

get_Omega_beta(G, [2*H + 4*c], 3; check_predictions=true)
# All predictions hold.