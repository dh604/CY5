import time
load("PT_CY5.sage")

# ---------- Our goal --------------------------------------------------
# We want to confirm the the claims made in Evidence 5.7 about the PT
# series of a threefold X in
# 
# Z1 = Tot_P1xP1 O(-1,-1)+O(-1,-1)+O(0,0)
# 
# and
# 
# Z2 = Tot_P1xP1 O(-2,0)+O(0,-2)+O(0,0)
#
# up to degree

dmax = 4

# and with precision

qprec = 7

# in the box-counting variable q. In particular, we want to verify
# Conjecture 5.6 which claims that the PT series of any threefold
# X in Z1 vanishes. As explained in Evidence 5.7 we should also
# expect a vanishing for Z2 once multicovering contributions of
# curves in class (1,0) and (0,1) are subtracted.



# ---------- The PT series ---------------------------------------------
# First we want to define a function which computes the PT series of
#
# X = Tot_P1xP1 O(a1,a2)
#
# in
#
# Z = Tot_P1xP1 O(a1,a2)+O(b1,b2)+O(c1,c2)
#
# for arbitrary a, b, c. First we set up the ring for our calculation:

setup_ring('q1,q2,q3,q4')

# and introduce an auxiliary variable

q5 = 1/(q1*q2*q3*q4)

# We evaluate the PT series via the K-theoretic vertex:

def PT_local_P1xP1(dd, LB_degs, qprec):
    # The tangent space of Z at the fixed points decomposes into the following
    # T representations. These are the square roots of the characters mentioned
    # in the paper since Henry's code for the vertex doubles all characters anyway.
    [a, b, c] = LB_degs
    
    wts = [
        [q1, q2, q3, q4, q5],
        [q2, 1/q1, q3*q1^(-a[0]), q4*q1^(-b[0]), q5*q1^(-c[0])],
        [1/q1, 1/q2, q3*q1^(-a[0])*q2^(-a[1]), q4*q1^(-b[0])*q2^(-b[1]), q5*q1^(-c[0])*q2^(-c[1])],
        [1/q2, q1, q3*q2^(-a[1]), q4*q2^(-b[1]), q5*q2^(-c[1])]
    ]
    
    # The normal bundle of the compact edges splits into line bundles
    # with degrees either of the following two:
    LB_degs_loc = [
        [0, a[0], b[0], c[0]],
        [0, a[1], b[1], c[1]]
    ]
    
    # Now we compute the PT invariant via the K-theoretic vertex
    return sum(V([mu[0], nu[0], []], wts[0], qprec) *
               V([nu[1], mu[0].conjugate(), []], wts[1], qprec) *
               V([mu[1].conjugate(), nu[1].conjugate(), []], wts[2], qprec) *
               V([nu[0].conjugate(), mu[1], []], wts[3], qprec) *
               E(mu[0], LB_degs_loc[0], wts[0]) *
               E(nu[1], LB_degs_loc[1], wts[1]) *
               E(mu[1].conjugate(), LB_degs_loc[0], wts[2]) *
               E(nu[0].conjugate(), LB_degs_loc[1], wts[3])
               for mu in partition_list(dd[0], 2)
               for nu in partition_list(dd[1], 2)
               )

# Since we will only require the Calabi-Yau speicalisation a + b + c = (-2,-2)
# we introduce the following function

PT_local_P1xP1_CY = lambda dd, aa, qprec: PT_local_P1xP1(dd, [[-aa[0],-aa[1]], [-aa[2],-aa[3]], [aa[0]+aa[2]-2,aa[1]+aa[3]-2]], qprec)


# ---------- Tot_P1xP1 O(-1,-1)+O(-1,-1)+O(0,0) ------------------------
# We want to test Conjecture 5.6 which predicts the vanishing of all PT
# invariants of any threefold in the above fivefold. We recall this goal
# here:

print("The PT invariants of the following threefolds X in Tot_P1xP1 "
      "O(a1,a2) + O(b1,b2) + O(c1,c2) all cojecturally vanish. We always "
      "choose X = Tot_P1xP1 O(a1,a2).\nWe test this conjecture up to {} "
      "orders beyond leading order in q in q and bi-degrees (d1,d2) with "
      "d1 + d2 <= {}:".format(qprec - 1, dmax))

# By obvious symmetries, the only two cases to check are
#
#  X = Tot_P1xP1 O(-1,-1)   and   X = Tot_P1xP1 O(0,0)
#
# which we check in the below loop:

for aa in [[1, 1, 1, 1], [0, 0, 1, 1]]:
    
    print("\n" + 62 * "-" + 
          "\nTot_P1xP1 O({},{})  in  Tot_P1xP1 O({},{}) + O({},{}) + O({},{})"
          "\n".format(-aa[0], -aa[1], -aa[0], -aa[1], -aa[2], -aa[3], aa[0]+aa[2]-2, aa[1]+aa[3]-2) +
          62 * "-")
    
    # Now we check the vanishing in all bi-degrees (d1,d2)
    # with d1 + d2 <= dmax:
    
    for d1 in range(dmax + 1):
        for d2 in range(dmax + 1 - d1):
            
            # We evaluate the PT series
            
            start = time.time()
            PT = PT_local_P1xP1_CY([d1, d2], aa, qprec)
            
            # it always vanishes except when d1=d2=0 where
            # by definition the series is 1. So the conjectural
            # expression to match is
            
            if d1 + d2 ==0:
                conj = 1
            else:
                conj = 0
            
            # Now we test the vanishing
            
            conj_holds = (PT - conj).is_zero()
            end = time.time()
            
            # and report the result
            
            print("(d1,d2) = ({},{})  -->  {}  "
                  "({:.2f}h)".format(d1, d2, conj_holds, (end - start) / 360.)
            )
            
            # If the conjecture is false we print the
            # discrepancy:
            
            if conj_holds == False:
                print("Difference PT series - conjectural formula =")
                print((PT, conj))

# ---------- Tot_P1xP1 O(-2,0)+O(0,-2)+O(0,0) ------------------------
# We repeat this for the above fivefold. However, now the conjecture is
# not that the PT series equally vanish but that they equate to
#
# Exp( Q1 * Omega1 + Q2 * Omega2 )
#
# for some suitable Omega1, Omega2. See Evidence 5.7. So first we compute
# the plethystic exponential.

print("\nStart setting up the conjectural formulas for the next PT series "
      "to analyse...")

# We compute Exp( ... ) in the ring

QR.<Q1,Q2> = PowerSeriesRing(wR.fraction_field(), default_prec = dmax + 1)

start = time.time()

# The argument of Exp is

Omega_P1xP1_20_02_00 = (
    - Q1 * ((q2^2*(-1+q1^2*q3^2)*q4^2*(1+q2^2*q4^2))/((-1+q4^2)*(-1+q2^4*q4^2)*(-1+q1^2*q2^2*q3^2*q4^2)))
    - Q2 * (q1^2*q3^2*(1+q1^2*q3^2)*(-1+q2^2*q4^2))/((-1+q3^2)*(-1+q1^4*q3^2)*(-1+q1^2*q2^2*q3^2*q4^2))
)

# Now we expand up to order dmax and extract the coefficients

PT_P1xP1_20_02_00_dict = Exp(Omega_P1xP1_20_02_00, dmax).coefficients()

# So far the PT series are presented as rational functions. With
# the following function we replace q4 with q*q4 and perform a
# Laurent expansion in q:

q_expand = Hom(wR.fraction_field(), qR)([q1, q2, q3, q4*q])

# By the symmetries of the fivefold we only need to check our
# prediction for the threefolds
#
#  X1 = Tot_P1xP1 O(-2,0)   and   X2 = Tot_P1xP1 O(0,0) .
#
# To get a prediction for the first threefold we just need to
# apply q_expand to PT_P1xP1_20_02_00_dict. To get a prediction
# for X2 we have to swap q3 and q5. The following function does
# exactly that and provides us with a prediction for the PT
# series of X1 and X2.

def get_PT_P1xP1_conj(d1, d2, aa):
    if aa == [2, 0, 0, 2]:    # if X1 then
        return q_expand(PT_P1xP1_20_02_00_dict[Q1^d1 * Q2^d2])
    elif aa == [0, 0, 0, 2]:  # if X2 then
        swap_q3_q5 = Hom(wR.fraction_field(), wR.fraction_field())([q1, q2, q5, q4])
        return q_expand(swap_q3_q5(PT_P1xP1_20_02_00_dict[Q1^d1 * Q2^d2]))

end = time.time()

print("Done! This took {:.2f}h.\nNow we check the formula up to "
      "{} orders beyond leading order in q and bi-degrees (d1,d2) "
      "with d1 + d2 <= {}:".format((end - start)/360., qprec - 1, dmax))

# Now that we have our prediction, let's test it:

for aa in [[2, 0, 0, 2], [0, 0, 0, 2]]:  # for X1 and X2:   
    
    print("\n" + 59 * "-" +
          "\nTot_P1xP1 O({},{})  in  Tot_P1xP1 O({},{}) + O({},{}) + "
          "O({},{})\n".format(-aa[0], -aa[1], -aa[0], -aa[1], -aa[2], -aa[3], aa[0]+aa[2]-2, aa[1]+aa[3]-2) +
          59 * "-")
    
    # Test the prediction in bi-degree (d1,d2) with
    # d1 + d2 <= dmax.
    
    for d1 in range(dmax + 1):
        for d2 in range(dmax + 1 - d1):
            
            start = time.time()
            
            # The PT series is
            
            PT = PT_local_P1xP1_CY([d1, d2], aa, qprec)
            
            # the conjectural formula is
            
            conj = get_PT_P1xP1_conj(d1, d2, aa)
            
            # Check if they match
            
            conj_holds = (PT - conj).is_zero()
            end = time.time()
            
            # and report the result.

            print("(d1,d2) = ({},{})  -->  {}  "
                  "({:.2f}h)".format(d1, d2, conj_holds, (end - start) / 360.))
            
            # If they mismatch, print the difference.
            
            if conj_holds == False:
                print("Difference PT series - conjectural formula =")
                print((PT, conj))