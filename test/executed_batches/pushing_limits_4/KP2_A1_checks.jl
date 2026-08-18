using Oscar, GKMtools, CY5

# Test conjecture for Tot(K_{P2}) x A_r

r = 1

Y = projective_space(GKM_graph, 2)

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")

# curve class of P2
H = curve_class(G, "1,1", "2,1")

gromov_witten(G, 3*H + b, 0, class_one(); g=0)
# 0

gromov_witten(G, 3*H + b, 0, class_one(); g=1)
# 0

gromov_witten(G, 3*H + b, 0, class_one(); g=2)
# 360*t4^2