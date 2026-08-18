using Oscar, GKMtools, GW_CY5

N = 3
m = -1

G = gkm_5d_gauge(N, m; equiCY=true)
F = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]
B = [curve_class(G, "v$i", "w$i") for i in 1:N]

gMax = 2

get_Omega_beta(G, [2*F[1] + 2*F[2]], gMax; check_predictions=true)
# All predicitons hold.

get_Omega_beta(G, [2*(B[2] + F[1] + F[2])], gMax; check_predictions=true)
# Killed.

get_Omega_beta(G, [3*F[1], 4*F[1]], gMax; check_predictions=true)