# Test conjecture for A_s x C^1 x A_r

s = 2
r = 3

As_times_C1 = Ar_times_C1(s)
G = X_times_Ar(As_times_C1, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of P2
B = curve_class(G, "1,1", "2,1")
C = curve_class(G, "2,1", "3,1")

# run tests (up to genus 2)

get_Omega_beta(G, [2*b], 2; check_predictions=true)
# zero for 2*b, nonzero for b

get_Omega_beta(G, [2*B], 2; check_predictions=true)
# zero for 2*b, nonzero for b

get_Omega_beta(G, [2*B + 2*b], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [B + 2*b], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [2*B + b], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [3*B + 2*b], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [2*b+2*c], 2; check_predictions=true)
# nonzero for b+c, zero for 2(b+c)

get_Omega_beta(G, [b+C], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [b+B+C], 2; check_predictions=true)
# all zero

get_Omega_beta(G, [4*c], 2; check_predictions=true)
# c is nonzero, 2c, 4c are zero