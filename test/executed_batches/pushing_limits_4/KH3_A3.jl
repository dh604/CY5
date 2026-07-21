# Test conjecture for Tot(K_{Hirzebruch_3}) x A_r

using Oscar, GKMtools, CY5

r = 3

Y = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 3))

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of P2
H1 = curve_class(G, "1,1", "2,1") # Chern number -1 in Y
H2 = curve_class(G, "1,1", "4,1") # Chern number 2 in Y (i.e. fiber)

# run tests (tested some up to genus 3, the others up to genus 2.)

get_Omega_beta(G, [4*b, 3*b], 3; check_predictions=true)
# holds for 2*b, 3*b, 4*b, fails for b.

get_Omega_beta(G, [4*H1], 2; check_predictions=true)
# fails for 4*H1, 2*H1, H1.

get_Omega_beta(G, [2*H2], 3; check_predictions=true)
# fails for H2, holds for 2*H2

get_Omega_beta(G, [2*(b+H1)], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*(b+c)], 3; check_predictions=true)
# holds for 2*(b+c), fails for b+c.

get_Omega_beta(G, [2*(H1+H2)], 2; check_predictions=true)
# fails for H1+H2, 2*(H1+H2).

get_Omega_beta(G, [2*b + H1], 3; check_predictions=true)
# All predictions hold.

##

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

get_Omega_beta(G, [H2 + 2*b], 3; check_predictions=true)

get_Omega_beta(G, [2*H2 + 2*b], 3; check_predictions=true)

get_Omega_beta(G, [H2 + 2*c], 3; check_predictions=true)