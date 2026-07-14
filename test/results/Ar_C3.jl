# File to test A_r x C^3.

using Oscar, GKMtools, CY5

r = 3
G = Ar_times_C3(r)
b = curve_class(G, Edge(1, 2))
c = curve_class(G, Edge(2, 3))
d = curve_class(G, Edge(3 ,4))

max_genus = 3

get_Omega_beta(G, [6*b], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [3*b + 3*c], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [5*b], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [4*b + 2*c], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*(b+c+d)], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [2*b+c+d], max_genus; check_predictions=true)
# All predictions hold.

get_Omega_beta(G, [b+2*c+d], max_genus; check_predictions=true)
# All predictions hold.

################# HAVE THE FOLLOWING FINISHED ? #################
#
#       (part of the pushing_limits_3/ batch)
#
#################################################################

get_Omega_beta(G, [4*b+4*c+4*d, 3*b+3*c+3*d], max_genus; check_predictions=true)

get_Omega_beta(G, [4*b+2*c+4*d], max_genus; check_predictions=true)

get_Omega_beta(G, [4*b+4*c+2*d], max_genus; check_predictions=true)