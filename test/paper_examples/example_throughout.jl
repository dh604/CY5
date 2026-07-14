##################################
# Example 1.2
##################################

using GKMtools, CY5
G = gkm_graph_of_example_1_1();
beta = curve_class(G, Edge(1, 2));
gMax = 1;
GW = get_GW_beta(G, [beta, 2*beta, 3*beta], gMax; show_bar=false);
GW[beta]
GW[2*beta]
GW[3*beta]

# julia> GW[beta]
# (1//12*t1*t4*t5*u^2 + 1//12*t1*t4*t6*u^2 + 1//12*t1*t5*t6*u^2 + t1 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t4*t6*u^2 + 1//12*t2*t5*t6*u^2 + t2 + 1//12*t3*t4*t5*u^2 + 1//12*t3*t4*t6*u^2 + 1//12*t3*t5*t6*u^2 + t3)//(t4*t5*t6*u^2)

# julia> GW[2*beta]
# (1//24*t1*t4*t5*u^2 + 1//24*t1*t4*t6*u^2 + 1//24*t1*t5*t6*u^2 + 1//8*t1 + 1//24*t2*t4*t5*u^2 + 1//24*t2*t4*t6*u^2 + 1//24*t2*t5*t6*u^2 + 1//8*t2 + 1//24*t3*t4*t5*u^2 + 1//24*t3*t4*t6*u^2 + 1//24*t3*t5*t6*u^2 + 1//8*t3)//(t4*t5*t6*u^2)

# julia> GW[3*beta]
# (1//36*t1*t4*t5*u^2 + 1//36*t1*t4*t6*u^2 + 1//36*t1*t5*t6*u^2 + 1//27*t1 + 1//36*t2*t4*t5*u^2 + 1//36*t2*t4*t6*u^2 + 1//36*t2*t5*t6*u^2 + 1//27*t2 + 1//36*t3*t4*t5*u^2 + 1//36*t3*t4*t6*u^2 + 1//36*t3*t5*t6*u^2 + 1//27*t3)//(t4*t5*t6*u^2)

##################################
# Example 1.4
##################################

using GKMtools, CY5
G = gkm_graph_of_example_1_1();
beta = curve_class(G, Edge(1, 2));
gMax = 1;
Omega = get_Omega_beta(G, [beta, 2*beta, 3*beta], gMax; show_bar=false);
Omega[beta]
Omega[2*beta]
Omega[3*beta]

# julia> Omega[beta]
# (1//12*t1*t4*t5*u^2 + 1//12*t1*t4*t6*u^2 + 1//12*t1*t5*t6*u^2 + t1 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t4*t6*u^2 + 1//12*t2*t5*t6*u^2 + t2 + 1//12*t3*t4*t5*u^2 + 1//12*t3*t4*t6*u^2 + 1//12*t3*t5*t6*u^2 + t3)//(t4*t5*t6*u^2)

# julia> Omega[2*beta]
# 0

# julia> Omega[3*beta]
# 0