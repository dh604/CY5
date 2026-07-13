# Test conjecture for Tot(K_{Hirzebruch_5}) x A_r

r = 3

Y = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5))

Ktot = total_space(wedge_product(cotangent_bd(Y), 2))
G = X_times_Ar(Ktot, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of P2
H1 = curve_class(G, "1,1", "2,1")
H2 = curve_class(G, "1,1", "4,1")

# run tests (tested some up to genus 3, the others up to genus 2.)

get_Omega_beta(G, [4*b], 3; check_predictions=true)
# nonzero for b, zero for 2b, 4b

get_Omega_beta(G, [4*H1], 2; check_predictions=true)
# nonzero for H1 and 2*H1

get_Omega_beta(G, [2*H2], 3; check_predictions=true)
# nonzero for H2, zero for 2*H2.

get_Omega_beta(G, [2*(b+H1)], 3; check_predictions=true)
# get zero for b+H1, get nonzero for 2(b+H1)

get_Omega_beta(G, [2*(b+c)], 3; check_predictions=true)
# get zero for 2(b+c), nonzero for b+c

get_Omega_beta(G, [2*(H1+H2)], 2; check_predictions=true)
# both nonzero

get_Omega_beta(G, [2*b + H1], 3; check_predictions=true)
# get zero