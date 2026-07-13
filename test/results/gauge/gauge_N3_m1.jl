using Oscar, GKMtools, CY5

N = 3
m = -1

G = gkm_5d_gauge(N, m; equiCY=true)
F = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]
B = [curve_class(G, "v$i", "w$i") for i in 1:N]

gMax = 2

################# HAVE THE FOLLOWING FINISHED ? #################
#
#       (part of the pushing_limits_2/ batch)
#
#################################################################

get_Omega_beta(G, [B[1]], gMax; check_predictions=true)

get_Omega_beta(G, [2*B[2]], gMax; check_predictions=true)

get_Omega_beta(G, [3*F[1], 4*F[1]], gMax; check_predictions=true)

get_Omega_beta(G, [2*F[2] + 2*B[3]], gMax; check_predictions=true)

get_Omega_beta(G, [4*B[3], 3*B[3]], gMax; check_predictions=true)

get_Omega_beta(G, [2*F[1] + 2*F[2]], gMax; check_predictions=true)

get_Omega_beta(G, [2*(B[2] + F[1] + F[2])], gMax; check_predictions=true)