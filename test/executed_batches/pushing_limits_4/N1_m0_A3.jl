using Oscar, GKMtools, CY5

# Test conjecture for X_{N,m} x A_r with N=1, m=0

r = 3

X = gkm_3d_gauge(1, 0)
G = X_times_Ar(X, r; equiCY = true)

# curve classes of A_r
b = curve_class(G, "v1,1", "v1,2")
c = curve_class(G, "v1,2", "v1,3")
d = curve_class(G, "v1,3", "v1,4")

# curve class of X_{N,m}
B = curve_class(G, "v1,1", "w1,1")


get_Omega_beta(G, [4*B], 3; check_predictions=true)
# holds for 2B, 4B, fails for B

get_Omega_beta(G, [B + b, 2*B + b, 3*B + b], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [B + c, 2*B + c, 3*B + c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [B + b+c, 2*B + b+c, 3*B + b+c], 3; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [B + 2*b, 2*B + 2*b, 3*B + 2*b], 3; check_predictions=true)
# All predictions hold.

# Added later, not in file:
get_Omega_beta(G, [B+b+c+d, 2*B+b+c+d, 2*(B+b+c)], 3; check_predictions=true)
# All predictions hold.

# Added later, still running:
get_Omega_beta(G, [3*(B+b+c)], 3; check_predictions=true)
# Killed.