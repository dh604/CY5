# Test conjecture for Tot(K_{P2}) x A_r

using Oscar, GKMtools, GW_CY5

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

# run tests

get_Omega_beta(G, [4*b], 2; check_predictions=true)
# Get nonzero stuff starting in genus zero for b, but zero for 2b and 4b. OK as moduli space non-compact.

get_Omega_beta(G, [3*b], 2; check_predictions=true)
# Get nonzero stuff starting in genus zero for b, but zero for 3b. OK as moduli space non-compact.

get_Omega_beta(G, [2*H], 2; check_predictions=true)
# Get nonzero stuff starting in genus zero. OK as moduli space non-compact.

get_Omega_beta(G, [3*H, 4*H], 2; check_predictions=true)
# Get nonzero stuff starting in genus zero for H, 2H, 3H, 4H. OK as moduli space non-compact.

get_Omega_beta(G, [2*(H+b)], 3; check_predictions=true)
# Get zero up to genus 3.

get_Omega_beta(G, [2*H+b], 3; check_predictions=true)
# Get zero up to genus 3.

get_Omega_beta(G, [H+2*b], 3; check_predictions=true)
# Get zero up to genus 3.

get_Omega_beta(G, [(H + b + c)], 3; check_predictions=true)
# Get zero up to genus 3.

get_Omega_beta(G, [H + b + 2*c], 3; check_predictions=true)
# Get zero up to genus 3.

get_Omega_beta(G, [2*H + 4*c], 3; check_predictions=true)
# Get zero up to genus 3.