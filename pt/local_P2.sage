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
#    Z1 = Tot_P2 O(-1) + O(-1) + O(-1)
#
# partly stated in equation (30) and the remaining indices listed in the
# folder
#
#    omega_data_folder
#
# specified above. Moreover, we will test Conjecture 5.6 claiming that
# the PT series of any threefold inside
#
#    Z2 = Tot_P2 O(-2) + O(-1) + O(0)
#
# vanishes.
#
#
# ---------- The PT series ---------------------------------------------
# For this we first want to define a function which computes the PT
# series of
#
# X = Tot_P2 O(a)
#
# in
#
# Z = Tot_P2 O(a) + O(b) + O(c)
#
# for arbitrary a, b, c. First we set up the ring for our calculation:

setup_ring('q0,q1,q2,q3,q4')

# and introduce an auxiliary variable

q5 = 1/(q0*q1*q2*q3*q4)

# The following function evaluates the degree d PT series of Z with
# line bundle degrees LBdegs = [a,b,c]. When a=b=c=-1 the decomposition
# of the tangent space of Z at a fixed point coincides with the
# convention used in Section 4.1.1.

def PT_local_P2(d, LBdegs, qprec):

    # The tangent space of Z at the fixed points decomposes into the following
    # T' representations. These are the square roots of the characters mentioned
    # in the paper since Henry's code for the vertex doubles all characters anyway.
    wts = [
        [q1/q0, q2/q0, q3*q0^(-LBdegs[0]), q4*q0^(-LBdegs[1]), q5*q0^(-LBdegs[2])],
        [q2/q1, q0/q1, q3*q1^(-LBdegs[0]), q4*q1^(-LBdegs[1]), q5*q1^(-LBdegs[2])],
        [q0/q2, q1/q2, q3*q2^(-LBdegs[0]), q4*q2^(-LBdegs[1]), q5*q2^(-LBdegs[2])]
    ]

    # The normal bundle of all three edges splits into line bundles of degrees
    LBdegs = [1] + LBdegs

    # Now we compute the PT invariant via the K-theoretic vertex
    return sum(
               prod(
                    V([mu[i], mu[i - 1 % 3].conjugate(), []], wts[i], qprec) *
                    E(mu[i], LBdegs, wts[i])
                    for i in range(3)
                    )
               for mu in partition_list(d, 3)
               )



# ---------- Tot_P2 O(-1) + O(-1) + O(-1) ------------------------------
# Note that by the symmetries of this fivefold we only need to analyse
# the PT series of one threefold, namely
#
#   X = Tot_P2 O(-1) .
#

print(28 * "-" + "\nTot_P2 O(-1) + O(-1) + O(-1)\n" + 28 * "-")

# We begin by reading the conjectural formulas for the membrane indices
# from the files.

start = time.time()

with open(omega_data_folder + 'Omega_P2_111_num.dat') as f:
    Omega_num = f.read().splitlines()

with open(omega_data_folder + 'Omega_P2_111_den.dat') as f:
    Omega_den = f.read().splitlines()

# We check conjectural formulae for membrane indices against the PT series
# up to degree

dmax = len(Omega_num)

# and precision

qprec = 8

# in the box-counting variable q. To get a prediction for the PT series we
# need to compute Exp ( sum_d Q^d Omega_ d ). This calculation is done in
# the ring

QR.<Q> = PowerSeriesRing(wR.fraction_field(), default_prec = dmax + 1)

# We denote the argument of the plethystic exponential by

Omega_P2_O111 = sum(QR(Omega_num[i]) / QR(Omega_den[i]) * Q^(i+1)
                    for i in range(dmax)) + O(Q^(dmax+1))
end = time.time()

print("Finished reading Omega. This took {:.2f}h.\nNow we compute the "
      "plethystic exponential".format((end - start)/360.))

# To match with the PT series we shift q4 -> q*q4 and perform
# a Laurent expansion in q via the following function

q_expand = Hom(wR.fraction_field(), qR)([q0, q1, q2, q3, q4*q])

# Let us now compute the conjectural expressions for the PT series
# and store them in a list.

start = time.time()
PT_P2_O111_conj = Exp(Omega_P2_O111, dmax).padded_list()
end = time.time()

print("Expansion of conjectural formula took {:.2f}h.\nNow we check "
      "the formula up to {} orders beyond leading order in q and degree "
      "d <= {}.\n".format((end - start)/360., qprec - 1, dmax))

# Now we can compare the PT series calculated by PT_local_P2 via the
# K-theoretic vertex with our conjectural formulae stored in PT_P2_O111_conj.

for d in range(1, dmax + 1):
    start = time.time()
    LBdegs = [-1, -1, -1]
    PT = PT_local_P2(d, LBdegs, qprec)
    conj = q_expand(PT_P2_O111_conj[d])
    conj_holds = (PT - conj).is_zero()
    end = time.time()
    print("Tot_P2 O({}) + O({}) + O({}), d = {}  -->  {}  "
          "({:.2f}h)".format(*LBdegs, d, conj_holds, (end - start) / 360.))


# ---------- Tot_P2 O(-2) + O(-1) + O(0) -------------------------------
# Next we probe Conjecture 5.6 which claims the vanishing of the PT series
# of any threefold inside the above fivefold. Note that there are three
# threefolds to check:
#
#   X1 = Tot_P2 O(-2),  X1 = Tot_P2 O(-1)  and  X3 = Tot_P2 O(0).
#
# but there is also the choice of which of the remaining two line bundles is
# L4 and which is L5. In other words, we need to choose a line bundle on
# which the torus C^*_q acts with weight +1. This leaves us with 3*2 = 6
# different configurations to check.
# We will check the vanishing of all three series up to degree

dmax = 4

# with precision

qprec = 7

# in q. We check all 6 configurations in the loop below:

for LBdegs in Permutations([-2, -1, 0]):
    print("\n" + 45 * "-" +
          "\nTot_P2 O({})  in  Tot_P2 O({}) + O({}) + O({})\n".format(LBdegs[0], *LBdegs) +
          45 * "-")

    print("We test the vanishing of the PT series up to {} orders beyond "
          "leading order in q and degree d <= {}:".format(qprec - 1, dmax))

    LBdegs = list(LBdegs)
    for d in range(1, dmax + 1):
        start = time.time()
        PT = PT_local_P2(d, LBdegs, qprec)
        conj = 0
        conj_holds = (PT - conj).is_zero()
        end = time.time()
        print("Tot_P2 O({}) + O({}) + O({}), d = {}  -->  {}  "
              "({:.2f}h)".format(*LBdegs, d, conj_holds, (end - start) / 360.))
