# run 1 implemented.
using Oscar, GKMtools, CY5

# File to test A_r x C^3.

r = 3
G = Ar_times_C3(r)
b = curve_class(G, Edge(1, 2))
c = curve_class(G, Edge(2, 3))
d = curve_class(G, Edge(3 ,4))

max_genus = 1

get_Omega_beta(G, [2*b], max_genus; check_predictions=true)

get_Omega_beta(G, [2*b+c+d], max_genus; check_predictions=true)