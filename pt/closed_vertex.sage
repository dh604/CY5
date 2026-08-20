load("PT_CY5.sage")
import time


# ---------- Our goal --------------------------------------------------
# We want to confirm the formula in Conj 2.11 in bi-degree (d1,d2,d3)
# with d1 + d2 + d3 <= dmax with

dmax = 5

# and with precision

qprec = 7

# in the box-counting variable q.


# ---------- Setup of the PT series ------------------------------------
# First we set up the ring for the computation of the PT series:

setup_ring('q1,q2,q3')

# We introduce an auxiliary variable

kappa = 1/(q1*q2*q3)

# We choose the 𝖳′ characters at the fixed points of the GKM graph as

chars = [
    [1/(q1^2 * kappa), 1/(q3^2 * kappa), 1/(q2^2 * kappa), 1, kappa],
    [q1^2 * kappa, q3^2, q2^2, 1, kappa],
    [q3^2 * kappa, q2^2, q1^2, 1, kappa],
    [q2^2 * kappa, q1^2, q3^2, 1, kappa]
]

# where we doubled all weights in order to avoid fractional exponents.
# We also need the following cyclic permutations of chars[0] as arguments
# for edge terms:

chars01 = [1/(q3^2 * kappa), 1/(q2^2 * kappa), 1/(q1^2 * kappa), 1, kappa]
chars02 = [1/(q2^2 * kappa), 1/(q1^2 * kappa), 1/(q3^2 * kappa), 1, kappa]

# The normal bundles of the three torus invariant P1s split into line
# bundles of the following degrees:

LBdegs = [-1, -1, 0, 0]

# Now we are prepared to define a function which evaluates the PT series
# in curve class (d1,d2,d3) with precision qprec in q:

def PT_clvert(d1, d2, d3, qprec):
    return sum(
        V([mu1, mu3, mu2], chars[0], qprec) *
        V([mu1.conjugate(), [], []], chars[1], qprec) *
        V([mu3.conjugate(), [], []], chars[2], qprec) *
        V([mu2.conjugate(), [], []], chars[3], qprec) *
        E(mu1, LBdegs, chars[0]) *
        E(mu3, LBdegs, chars01) *
        E(mu2, LBdegs, chars02)
        for mu1 in Partitions(d1)
        for mu2 in Partitions(d2)
        for mu3 in Partitions(d3)
    )



# ---------- The conjectural formula -----------------------------------
# The conjectural formula is presented as a plethystic exponential. Since
# we would only like to verify the formula in multi-degree (d1,d2,d3)
# with d1 + d2 + d3 <= dmax, we evalute the plethystic exponential up to
# this order and extract conjectural formulae for each individual PT series.
# First we set up the ring in which the plethystic exponential will take
# values in:

wR2.<qq1,qq2,qq3,qq4> = LaurentPolynomialRing(QQ)
QR2.<Q1,Q2,Q3> = PowerSeriesRing(wR2.fraction_field(), default_prec = dmax + 1)

# We choose to compute the plethystic exponential in this ring first and
# only later perform the Laurent expansion in q4 = q in order to increase
# performance.

# We introduce auxiliary variables

kkappa = 1/(qq1*qq2*qq3)
qq5 = kkappa / qq4

# and the function [t] = t - t^(-1). 

def br(t):
    return t - t^(-1)

# The following function evaluates the plethystic exponential which
# according to Conjecture 2.11 equates to the PT series of the closed vertex:

def PT_clvert_conj_series(dmax):
    Omega100 = 1 / (br(qq4) * br(qq5))
    Omega110 = - br(qq3^2*kkappa) / (br(qq3^2) * br(qq4) * br(qq5))
    Omega101 = - br(qq2^2*kkappa) / (br(qq2^2) * br(qq4) * br(qq5))
    Omega011 = - br(qq1^2*kkappa) / (br(qq1^2) * br(qq4) * br(qq5))
    Omega111 = 1 / (br(qq4) * br(qq5))
    
    return Exp(
               Omega100 * (Q1 + Q2 + Q3) +
               Omega110 * Q1 * Q2 +
               Omega101 * Q1 * Q3 +
               Omega011 * Q2 * Q3 +
               Omega111 * Q1 * Q2 * Q3
    , dmax)

# Let us evaluate the conjectural expression for the PT series and record
# the formulae in a dictionary.

start = time.time()  
PT_clvert_conj_dict = PT_clvert_conj_series(dmax).coefficients()
end = time.time()

print(
    "Expansion of the plethystic exponential relating membrane indices "
    "and the PT series took {:.2f}h.\n"
    "Now we check the formulae against PT series in multi-degree "
    "(d1, d2, d3) with d1 + d2 + d3 <= {} for {} orders beyond the "
    "leading order in q.\n".format((end - start)/360., dmax, qprec - 1)
)

# The conjectural formulae for the PT series are still elements of the
# ring wR2. We perform a Laurent expansion in q4=q up order qprec using

def q_expand(f, qprec):
    f_expanded = Hom(wR2,qR)([q1, q2, q3, q]).extend_to_fraction_field()(f)
    return f_expanded + O(q^(f_expanded.valuation() + 2 * qprec))


# ---------- Comparison/Check of Conjecture 2.11 -----------------------
# Now we can systematically check that the PT series computed via PT_clvert
# indeed equate to the conjectural formulae stored in PT_clvert_conj_dict:

for d in range(1, dmax+1):
    for part in Partitions(d, max_length=3):
        dList = part + [0] * (3 - len(part))
        
        start = time.time()
        PT = PT_clvert(*dList, qprec)
        conj = q_expand(PT_clvert_conj_dict[Q1^dList[0] * Q2^dList[1] * Q3^dList[2]], qprec)
        conj_holds = (PT - conj).is_zero()
        end = time.time()
        
        print("[d1, d2, d3] = {}  -->  {}  ({:.2f}h)".format(dList, conj_holds, (end - start) / 360.))