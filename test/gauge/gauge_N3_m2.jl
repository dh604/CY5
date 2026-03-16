N = 3
m = 2

G = gkm_5d_gauge(N, m; equiCY=true)
Qt = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]
Q = [curve_class(G, "v$i", "w$i") for i in 1:N]

gMax = 2

get_Omega_beta(G, [Q[1]], gMax)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (3, 1, 1) => (1399//1152*t4^4*u^4 + 21//4*t4^3*t5*u^4 + 4661//576*t4^2*t5^2*u^4 + 115//24*t4^2*u^2 + 21//4*t4*t5^3*u^4 + 10*t4*t5*u^2 + 1399//1152*t5^4*u^4 + 115//24*t5^2*u^2 + 5)//(t4*t5*u^2)

get_Omega_beta(G, [2*Q[2]], gMax)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (0, 1, 1) => (29//640*t4^4*u^4 + 1//4*t4^3*t5*u^4 + 27//64*t4^2*t5^2*u^4 + 7//8*t4^2*u^2 + 1//4*t4*t5^3*u^4 + 2*t4*t5*u^2 + 29//640*t5^4*u^4 + 7//8*t5^2*u^2 + 3)//(t4*t5*u^2)
#   (0, 2, 2) => (-133//40*t4^4*u^4 - 14*t4^3*t5*u^4 - 171//8*t4^2*t5^2*u^4 - 17//2*t4^2*u^2 - 14*t4*t5^3*u^4 - 35//2*t4*t5*u^2 - 133//40*t5^4*u^4 - 17//2*t5^2*u^2 - 6)//(t4*t5*u^2)

get_Omega_beta(G, [3*Q[3], 2*Q[3]], 2)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 3 entries:
#   (0, 2, 0) => 0
#   (0, 1, 0) => (7//5760*t4^4*u^4 + 1//576*t4^2*t5^2*u^4 - 1//24*t4^2*u^2 + 7//5760*t5^4*u^4 - 1//24*t5^2*u^2 + 1)//(t4*t5*u^2)
#   (0, 3, 0) => 0

get_Omega_beta(G, [3*Qt[1], 2*Qt[2], 2*Qt[2], 2*(Qt[1]+Qt[2])], 2)
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 6 entries:
#   (2, 0, 2) => 0
#   (3, 0, 0) => 0
#   (0, 0, 1) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)
#   (1, 0, 1) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)
#   (0, 0, 2) => 0
#   (1, 0, 0) => (1//360*t4^4*u^4 - 1//72*t4^2*t5^2*u^4 - 1//6*t4^2*u^2 - 1//2*t4*t5*u^2 + 1//360*t5^4*u^4 - 1//6*t5^2*u^2 - 2)//(t4*t5*u^2)