using Oscar, GKMtools, GW_CY5

N = 2
m = -1

G = gkm_5d_gauge(N, m; equiCY=true)
F = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]
B = [curve_class(G, "v$i", "w$i") for i in 1:N]

gMax = 2 # All predictions hold.

get_Omega_beta(G, [B[1]], gMax; check_predictions=true)

get_Omega_beta(G, [4*B[2], 3*B[2]], gMax; check_predictions=true)

get_Omega_beta(G, [3*F[1], 4*F[1]], gMax; check_predictions=true)

get_Omega_beta(G, [2*F[1] + 2*B[2]], gMax; check_predictions=true)