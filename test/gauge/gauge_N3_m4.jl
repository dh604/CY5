# (N, m) = (3, 4), outside the range of the conjecture.
# Computed values show dependence other than t4 and t5.

N = 3
m = 4

G = gkm_5d_gauge(N, m; equiCY=true)
Qt = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]
Q = [curve_class(G, "v$i", "w$i") for i in 1:N]

gMax = 2

# get_Omega_beta(G, [Q[1]], gMax)
# about 1.5h even for genus zero, since m >> 0 causes more ways of expressing Q[1].

get_Omega_beta(G, [Q[2]], gMax)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (0, 1, 3) => (1399//1152*t4^4*u^4 + 21//4*t4^3*t5*u^4 + 4661//576*t4^2*t5^2*u^4 + 115//24*t4^2*u^2 + 21//4*t4*t5^3*u^4 + 10*t4*t5*u^2 + 1399//1152*t5^4*u^4 + 115//24*t5^2*u^2 + 5)//(t4*t5*u^2)

get_Omega_beta(G, [3*Q[3], 2*Q[3]], 2)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 3 entries:
#   (0, 2, 0) => (-13//360*t1^8*t4*u^4 - 13//360*t1^8*t5*u^4 - 91//90*t1^7*t3*t4*u^4 - 91//90*t1^7*t3*t5*u^4 - 14//15*t1^7*t4^2*u^4 - 28//15*t1^7*t4*t5*u^4 - 14//15*t1^7*t5^2*u^4 - 17837//1440*t1^6*t3^2*t4*u^4 - 17837//1…
#   (0, 1, 0) => (1//360*t1^5*t4*u^4 + 1//360*t1^5*t5*u^4 + 7//144*t1^4*t3*t4*u^4 + 7//144*t1^4*t3*t5*u^4 + 11//240*t1^4*t4^2*u^4 + 11//120*t1^4*t4*t5*u^4 + 11//240*t1^4*t5^2*u^4 + 41//120*t1^3*t3^2*t4*u^4 + 41//120*t1^3…
#   (0, 3, 0) => (3//10*t1^12*t4*u^4 + 3//10*t1^12*t5*u^4 + 63//5*t1^11*t3*t4*u^4 + 63//5*t1^11*t3*t5*u^4 + 34//3*t1^11*t4^2*u^4 + 68//3*t1^11*t4*t5*u^4 + 34//3*t1^11*t5^2*u^4 + 9697//40*t1^10*t3^2*t4*u^4 + 9697//40*t1^1…

get_Omega_beta(G, [3*Qt[1], 2*Qt[2], 2*Qt[2], 2*(Qt[1]+Qt[2])], 2)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 6 entries:
#   (2, 0, 2) => 0
#   (3, 0, 0) => 0
#   (0, 0, 1) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)
#   (1, 0, 1) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)
#   (0, 0, 2) => 0
#   (1, 0, 0) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)