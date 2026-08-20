import time
load("PT_CY5.sage")


# !!! Important !!!
# You possibly need to change the below string to the path where all
# predictions for membrane indices are stored in case you modified the
# path structure:

omega_data_folder = '../test/omega_data/'



# ---------- Our goal --------------------------------------------------
# We want to verify the formulae for membrane indices of
#
#    Z = Tot_P3 O(-a) + O(a-4)
#
# partly stated in Section 4.2 and the remaining indices listed in the
# folder
#
#    omega_data_folder
#
# specified above. We will check the formulae with precision

qprec = 8

# in q.


# ---------- The PT series ---------------------------------------------
# For this we first want to define a function which computes the PT
# series of
#
#    X = P3
#
# in
#
#    Z = Tot_P3 O(-a) + O(a-4)
#
# for arbitrary a. First we set up the ring for our calculation:

setup_ring('q0,q1,q2,q3,q4')

# and introduce an auxiliary variable

q5 = 1/(q0*q1*q2*q3*q4)

# The degree d PT series of P3 in Z with precision qprec in q is
# then evaluated by the following function. We adopted the weight
# conventions from Section 4.2.

def PT_local_P3(d, a, qprec):
    # The tangent space of Z at the fixed points decomposes into the following
    # T representations. These are the square roots of the characters mentioned
    # in the paper since Henry's code for the vertex doubles all characters anyway.
    wts = [
        [q1/q0, q2/q0, q3/q0, q4*q0^a, q5*q0^(4-a)],
        [q2/q1, q0/q1, q3/q1, q4*q1^a, q5*q1^(4-a)],
        [q0/q2, q1/q2, q3/q2, q4*q2^a, q5*q2^(4-a)],
        [q0/q3, q2/q3, q1/q3, q4*q3^a, q5*q3^(4-a)]
    ]
    # Three of the edge terms will require a permutation of the characters as
    # arguments:
    wts_E = [
        [q3/q0, q1/q0, q2/q0, q4*q0^a, q5*q0^(4-a)],
        [q3/q1, q2/q1, q0/q1, q4*q1^a, q5*q1^(4-a)],
        [q3/q2, q0/q2, q1/q2, q4*q2^a, q5*q2^(4-a)]
    ]

    # The normal bundle of all three edges splits into line bundles of degrees
    LBdegs = [1, 1, -a, a - 4]

    # Now we compute the PT invariant via the K-theoretic vertex
    return sum(V([mu[0], mu[2].conjugate(), mu[3]], wts[0], qprec) *
               V([mu[1], mu[0].conjugate(), mu[4]], wts[1], qprec) *
               V([mu[2], mu[1].conjugate(), mu[5]], wts[2], qprec) *
               V([mu[3].conjugate(), mu[5].conjugate(), mu[4].conjugate()], wts[3], qprec) *
               E(mu[0], LBdegs, wts[0]) *
               E(mu[1], LBdegs, wts[1]) *
               E(mu[2], LBdegs, wts[2]) *
               E(mu[3], LBdegs, wts_E[0]) *
               E(mu[4], LBdegs, wts_E[1]) *
               E(mu[5], LBdegs, wts_E[2])
               for mu in partition_list(d, 6)
               )


# ---------- Membrane indices ------------------------------------------
# Next we want a function which computes the plethystic exponential of the
# generating series of membrane indices in order to get conjectural formulae
# for the PT series. The plethystic exponential takes values in

QR.<Q> = PowerSeriesRing(wR.fraction_field(), default_prec = 4)

# The following function loads the generating series of membrane indices for
# a given degree a from our files.

def get_Omega_series(a):

    # We only store formulae for membrane indices for a<=2.
    # When a>2, one can just switch the role of L4 and L5
    # by swapping q4 <-> q5. This is done by q4_trafo defined
    # below.

    if a<=2:
        LBdegs = [a, 4 - a]
        q4_trafo = End(wR).identity()
    else:
        LBdegs = [4 - a, a]
        q4_trafo = Hom(wR, wR)([q0, q1, q2, q3, q5])

    # Now we read the numerator and denominator of the membrane
    # indices

    with open('./omega_data/Omega_P3_{}{}_num.dat'.format(*LBdegs)) as f:
        Omega_num = f.read().splitlines()

    with open('./omega_data/Omega_P3_{}{}_den.dat'.format(*LBdegs)) as f:
        Omega_den = f.read().splitlines()

    dmax = len(Omega_num)

    # For one line bundle degree a, we present the formulae for membrane
    # indices in terms of elementary symmetric polynomials. The function
    # below transforms these into polynomials in the qi.

    ewR.<e1,e2,e3,e4> = PolynomialRing(wR)
    expand_e = Hom(ewR, wR)([q0+q1+q2+q3, q0*q1+q0*q2+q1*q2+q0*q3+q1*q3+q2*q3, q0*q1*q2+q0*q1*q3+q0*q2*q3+q1*q2*q3, q0*q1*q2*q3])

    # Now we form sum_d Q^d Omega_d and return the result.

    start = time.time()
    Omega_local_P3 = sum(q4_trafo(expand_e(ewR(Omega_num[i]))) / q4_trafo(wR(Omega_den[i])) * Q^(i+1) for i in range(dmax)) + O(Q^(dmax+1))
    end = time.time()

    print("Finished reading membrane indices of Tot_P3 O({}) + O({}). This "
          "took {:.2f}h.".format(-a, a-4, (end - start)/360.))

    return Omega_local_P3


# We will require the function below in order to compute the Laurent
# expansion of the membrane indices in q.

q_expand = Hom(wR.fraction_field(), qR)([q0, q1, q2, q3, q4*q])

# Now we test our conjectural formulae for a in {0,1,2,3,4}:

for a in range(5):
    print("\n" + 20 * "-" +
          "\nTot_P3 O({}) + O({})\n".format(-a,-4 + a) +
          20 * "-")

    # Set up the generating series of membrane indices

    Omega_local_P3 = get_Omega_series(a)

    # The series goes up to curve degree

    dmax = Omega_local_P3.prec() - 1

    # Now compute the plethystic exponential:

    print("Now we compute the plethystic exponential up to "
          "degree {}.".format(dmax))

    start = time.time()
    PT_local_P3_conj = Exp(Omega_local_P3, dmax).padded_list()
    end = time.time()

    print("\nExpansion of conjectural formula took {:.2f}h."
          "\nNow we check the formula up to {} orders beyond leading "
          "order in q.\n".format((end - start)/360., qprec - 1))

    # Now compare the PT series computed via the vertex with
    # our prediction.

    for d in range(1, dmax + 1):

        start = time.time()
        PT = PT_local_P3(d, a, qprec)
        conj = q_expand(PT_local_P3_conj[d])
        conj_holds = (PT - conj).is_zero()
        end = time.time()

        print("Tot_P3 O({}) + O({}), d = {}  -->  {}  "
              "({:.2f}h)".format(-a, -4+a, d, conj_holds, (end - start) / 360.))

        if conj_holds == False:
            print("Difference PT series - conjectural formula =")
            print(PT - conj)
