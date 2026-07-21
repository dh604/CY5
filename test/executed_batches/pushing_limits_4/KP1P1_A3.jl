using Oscar, GKMtools, CY5

# Test conjecture for Tot(K_{P1 x P1}) x A_r

r = 3

P1 = projective_space(GKM_graph, 1)
Y = P1 * P1

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1,1", "1,1,2")
c = curve_class(G, "1,1,2", "1,1,3")
d = curve_class(G, "1,1,3", "1,1,4")

# curve class of P1 x P1
H1 = curve_class(G, "1,1,1", "2,1,1")
H2 = curve_class(G, "1,1,1", "1,2,1")


get_Omega_beta(G, [2*H1 + 4*c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*(H1+b)], 3; check_predictions=true)

get_Omega_beta(G, [2*H1+b], 3; check_predictions=true)

get_Omega_beta(G, [3*H1+b], 3; check_predictions=true)

get_Omega_beta(G, [H1+2*b], 3; check_predictions=true)

get_Omega_beta(G, [(H1 + b + c)], 3; check_predictions=true)

get_Omega_beta(G, [H1 + b + 2*c], 3; check_predictions=true)

get_Omega_beta(G, [2*H1+H2+b], 3; check_predictions=true)

get_Omega_beta(G, [2*H1+2*H2+b], 3; check_predictions=true)

get_Omega_beta(G, [2*H1 + 2*H2 + 2*b], 3; check_predictions=true)

