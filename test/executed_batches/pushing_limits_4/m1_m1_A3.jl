using Oscar, GKMtools, CY5

# Test conjecture for a O(-1)+O(-1) chain of length s x A_r

s = 2
r = 3

X = minus_one_minus_one_chain_3d(s)
G = X_times_Ar(X, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "1,1", "1,2")
c = curve_class(G, "1,2", "1,3")
d = curve_class(G, "1,3", "1,4")

# curve class of X
B = curve_class(G, "1,1", "2,1")
C = curve_class(G, "2,1", "3,1")


get_Omega_beta(G, [4*B, 2*(B+C), 3*b], 3; check_predictions=true)

get_Omega_beta(G, [B + b, 2*B + b, 3*B + b], 3; check_predictions=true)

get_Omega_beta(G, [B+C + b, 2*B + C + b, 3*B+C + b], 3; check_predictions=true)

get_Omega_beta(G, [2*B + 2*C + b, 3*B+2*C + b], 3; check_predictions=true)

get_Omega_beta(G, [B + b+c, 2*B + b+c, 3*B + b+c], 3; check_predictions=true)

get_Omega_beta(G, [B+C + b+c, 2*B+C + b+c, 3*B+C + b+c], 3; check_predictions=true)

get_Omega_beta(G, [B+C + 2*b, 2*B+C + 2*b, 3*B+C + 2*b], 3; check_predictions=true)
