using Oscar, GKMtools, GW_CY5

# Test conjecture for a O(-1)+O(-1) chain of length s x A_r in nonproper beta.

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


get_Omega_beta(G, [4*B], 3; check_predictions=true)
# fails for B, holds for 2*B, 4*B

get_Omega_beta(G, [2*(B+C)], 3; check_predictions=true)
# fails for B+C, holds for 2*(B+C)

get_Omega_beta(G, [2*b], 3; check_predictions=true)
# fails for b, holds for 2*b

get_Omega_beta(G, [2*(b+c)], 3; check_predictions=true)
# fails for b+c, holds for 2*(b+c)

# Not in file, run manually:
get_Omega_beta(G, [3*(b+c), 3*(B+C), 3*B, 3*b], 3; check_predictions=true)
# holds for 3(b+c), 3(B+C), 3B, 3b;
# fails for (b+c), (B+C), B, b.