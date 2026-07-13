# TODO: run with corrected prediction for 4*b.

G = gkm_5d_P2_111()
b = curve_class(G, Edge(1, 2))

get_Omega_beta(G, [b, 2*b, 3*b, 4*b], 2; check_predictions=true)

# Prediction holds for (2)
# Prediction fails for (4):
#   prediction = 679//16384*t1^2*u^2 + 10129//32768*t1*t2*u^2 + 10129//32768*t1*t3*u^2 + 5343//16384*t1*t4*u^2 + 3//8*t1*t5*u^2 + 679//16384*t2^2*u^2 + 10129//32768*t2*t3*u^2 + 5343//16384*t2*t4*u^2 + 3//8*t2*t5*u^2 + 679//16384*t3^2*u^2 + 5343//16384*t3*t4*u^2 + 3//8*t3*t5*u^2 + 9885//32768*t4^2*u^2 + 3//8*t4*t5*u^2 + 3//8*t5^2*u^2 + 32583//32768
# Prediction holds for (1)
# Prediction holds for (3)

# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 4 entries:
#   (2) => -1//8*t1^2*u^2 + 7//32*t1*t2*u^2 + 7//32*t1*t3*u^2 + 3//32*t1*t4*u^2 + 3//32*t1*t5*u^2 - 1//8*t2^2*u^2 + 7//32*t2*t3*u^2 + 3//32*t2*t4*u^2 + 3//32*t2*t5*u^2 - 1//8*t3^2*u^2 + 3//32*t3*t4*u^2 + 3//32*t3*t5*u^2 + 3//32*t4^2*u^2…
#   (4) => 3//8*t1*t2*u^2 + 3//8*t1*t3*u^2 + 3//8*t1*t4*u^2 + 3//8*t1*t5*u^2 + 3//8*t2*t3*u^2 + 3//8*t2*t4*u^2 + 3//8*t2*t5*u^2 + 3//8*t3*t4*u^2 + 3//8*t3*t5*u^2 + 3//8*t4^2*u^2 + 3//8*t4*t5*u^2 + 3//8*t5^2*u^2
#   (1) => 1//32*t1^2*u^2 - 1//32*t1*t2*u^2 - 1//32*t1*t3*u^2 + 1//32*t2^2*u^2 - 1//32*t2*t3*u^2 + 1//32*t3^2*u^2 - 1//8
#   (3) => -9//32*t1*t2*u^2 - 9//32*t1*t3*u^2 - 9//32*t1*t4*u^2 - 9//32*t1*t5*u^2 - 9//32*t2*t3*u^2 - 9//32*t2*t4*u^2 - 9//32*t2*t5*u^2 - 9//32*t3*t4*u^2 - 9//32*t3*t5*u^2 - 9//32*t4^2*u^2 - 9//32*t4*t5*u^2 - 9//32*t5^2*u^2


get_Omega_beta(G, [b, 2*b, 3*b, 4*b], 1; check_predictions=true)

# Prediction holds for (2)
# Prediction fails for (4):
#   prediction = 32583//32768
# Prediction holds for (1)
# Prediction holds for (3)

# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 4 entries:
#   (2) => 1//8
#   (4) => 0
#   (1) => -1//8
#   (3) => 0

get_Omega_beta(G, [b, 2*b], 4; check_predictions=true)

# Prediction holds for (2)
# Prediction holds for (1)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (2) => -17//360*t1^6*u^6 + 1739//7680*t1^5*t2*u^6 + 1739//7680*t1^5*t3*u^6 + 217//2560*t1^5*t4*u^6 + 217//2560*t1^5*t5*u^6 - 763//1536*t1^4*t2^2*u^6 - 701//1536*t1^4*t2*t3*u^6 - 61//384*t1^4*t2*t4*u^6 - 61//384*t1^4*t2*t5*u^6 - 763/…
#   (1) => 17//23040*t1^6*u^6 - 17//7680*t1^5*t2*u^6 - 17//7680*t1^5*t3*u^6 + 7//1536*t1^4*t2^2*u^6 + 1//512*t1^4*t2*t3*u^6 + 7//1536*t1^4*t3^2*u^6 - 1//192*t1^4*u^4 - 25//4608*t1^3*t2^3*u^6 - 1//512*t1^3*t2^2*t3*u^6 - 1//512*t1^3*t2*t3…


get_Omega_beta(G, [b, 3*b], 3; check_predictions=true)

# Prediction holds for (1)
# Prediction holds for (3)
# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 2 entries:
#   (1) => -1//192*t1^4*u^4 + 1//96*t1^3*t2*u^4 + 1//96*t1^3*t3*u^4 - 1//64*t1^2*t2^2*u^4 - 1//64*t1^2*t3^2*u^4 + 1//32*t1^2*u^2 + 1//96*t1*t2^3*u^4 - 1//32*t1*t2*u^2 + 1//96*t1*t3^3*u^4 - 1//32*t1*t3*u^2 - 1//192*t2^4*u^4 + 1//96*t2^3*…
#   (3) => 3//4*t1^3*t2*u^4 + 3//4*t1^3*t3*u^4 + 3//4*t1^3*t4*u^4 + 3//4*t1^3*t5*u^4 - 89//64*t1^2*t2^2*u^4 - 65//32*t1^2*t2*t3*u^4 - 41//32*t1^2*t2*t4*u^4 - 41//32*t1^2*t2*t5*u^4 - 89//64*t1^2*t3^2*u^4 - 41//32*t1^2*t3*t4*u^4 - 41//32*…

get_Omega_beta(G, [b], 5; check_predictions=true)

# All predictions hold.
# Dict{AbstractAlgebra.FPModuleElem{ZZRingElem}, Any} with 1 entry:
#   (1) => -31//322560*t1^8*u^8 + 31//80640*t1^7*t2*u^8 + 31//80640*t1^7*t3*u^8 - 47//46080*t1^6*t2^2*u^8 - 1//1536*t1^6*t2*t3*u^8 - 47//46080*t1^6*t3^2*u^8 + 17//23040*t1^6*u^6 + 79//46080*t1^5*t2^3*u^8 + 1//1024*t1^5*t2^2*t3*u^8 + 1//…