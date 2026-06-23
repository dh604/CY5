#TODO: Human check this file.


# ==============================================================================
# Gauge-theory prediction of Omega_beta for the CY5 space Z_{N,m} = X_{N,m} x A^2.
#
# This file implements, from scratch, the partition function
#
#     Z_{N,m}(q4, q5; Q)
#
# of the *revised* specification in `GAUGE_latex_specs.tex` (the GW-Gauge
# correspondence). It supersedes the earlier `gauge_prediction.jl` and
# `gauge_prediction_2.jl`, which followed an earlier, mistaken specification.
# It is intentionally independent of those two files; it only relies on generic,
# non-gauge helpers of the package (`cc_mobius`, `downward_close_ccs`,
# `curve_class`) and on Oscar.
#
# ------------------------------------------------------------------------------
# Mathematical outline (everything is taken from the tex specification)
# ------------------------------------------------------------------------------
# The conjecture states that the disconnected equivariant GW partition function of
# Z_{N,m} equals  Z_{N,m}(e^{u eps_4}, e^{u eps_5}; Q)  as a series in the effective
# curve classes Q^beta of X_{N,m}.  We therefore
#
#   1. build the truncated series Z_{N,m} as a dictionary  {H_2-class -> coeff},
#   2. take its (ordinary) logarithm to obtain the *connected* free energy,
#   3. substitute  r4 = e^{u t4 / 2},  r5 = e^{u t5 / 2}  (so q4 = e^{u t4},
#      q5 = e^{u t5}) and expand in u,
#   4. apply the package's canonical BPS/Moebius inversion `cc_mobius` to read off
#      Omega_beta.
#
# Conventions.
#   * Half-integer powers of q4, q5 occur, so we work over
#         K = Frac(QQ[r4, r5]),     q4 = r4^2,   q5 = r5^2,
#     and every q-monomial q4^a q5^b is stored as r4^{2a} r5^{2b}.
#   * Curve classes are first tracked in nonnegative *monoid coordinates* -- the
#     multiplicities (d_1,...,d_N,e_1,...,e_{N-1}) of the generators
#         [B_a] = curve_class(G, "v$a", "w$a")          (sections,   a = 1..N)
#         [F_c] = curve_class(G, "v$c", "v$(c+1)")       (fibres,     c = 1..N-1)
#     -- because these are always >= 0 and bounded, so the componentwise
#     truncation is unambiguous (the H_2 coordinates of an effective class can be
#     negative).  The connected free energy is computed in monoid coordinates and
#     only then pushed forward to H_2.  The monoid -> H_2 pushforward imposes all
#     relations in H_2^eff; performing it after the (ordinary) logarithm is
#     equivalent to performing it before, because the pushforward is a ring
#     homomorphism and hence commutes with log.
# ==============================================================================


# ------------------------------------------------------------------------------
# Partition utilities
# ------------------------------------------------------------------------------

# Transpose (conjugate) of a partition given as a weakly decreasing Vector{Int}.
function _gp3_transpose(mu::Vector{Int})
  isempty(mu) && return Int[]
  cols = zeros(Int, mu[1])
  for x in mu, i in 1:x
    cols[i] += 1
  end
  return cols
end

# |mu| = sum of the parts.
_gp3_abs(mu::Vector{Int}) = sum(mu; init = 0)

# ||mu||^2 = sum of the squares of the parts (the symbol `sqabs` of the tex).
_gp3_sqnorm(mu::Vector{Int}) = sum(x^2 for x in mu; init = 0)

# All partitions of n, listed as weakly decreasing Vector{Int}.
function _gp3_partitions_into!(res, prefix, remaining, maxpart)
  if remaining == 0
    push!(res, copy(prefix))
    return
  end
  for p in min(remaining, maxpart):-1:1
    push!(prefix, p)
    _gp3_partitions_into!(res, prefix, remaining - p, p)
    pop!(prefix)
  end
end

# All partitions of 0, 1, ..., nmax (the empty partition is included as `Int[]`).
function _gp3_partitions_upto(nmax::Int)
  res = Vector{Vector{Int}}()
  for n in 0:nmax
    _gp3_partitions_into!(res, Int[], n, n)
  end
  return res
end


# ------------------------------------------------------------------------------
# Monomials in K = Frac(QQ[r4, r5])
# ------------------------------------------------------------------------------

# r4^e4 * r5^e5 as an element of K (exponents may be negative).
_gp3_mono(K, r4, r5, e4::Int, e5::Int) = K(r4)^e4 * K(r5)^e5


# ------------------------------------------------------------------------------
# The building block  \tilde Z_mu(A, B)
#
#     \tilde Z_mu(A, B) = prod_{(i,j) in mu} ( 1 - A^{mu^t_j - i + 1} B^{j - mu_i} )^{-1}
#
# `Aexp`, `Bexp` are the (r4, r5)-exponent pairs of A and B, e.g. A = q4 = r4^2 is
# encoded as Aexp = (2, 0).
# ------------------------------------------------------------------------------
function _gp3_tildeZ(mu::Vector{Int}, Aexp::Tuple{Int,Int}, Bexp::Tuple{Int,Int},
                    r4, r5, K)
  res = K(1)
  mut = _gp3_transpose(mu)
  for i in 1:length(mu), j in 1:mu[i]
    ea = mut[j] - i + 1        # exponent of A
    eb = j - mu[i]             # exponent of B
    e4 = ea * Aexp[1] + eb * Bexp[1]
    e5 = ea * Aexp[2] + eb * Bexp[2]
    res = res // (K(1) - _gp3_mono(K, r4, r5, e4, e5))
  end
  return res
end


# ------------------------------------------------------------------------------
# One-leg weight L_a for the a-th node (1-indexed), as an element of K.
#
# From the tex, the a-th factor of Z_{N,m} is
#     (-Q^{[B_a]})^{|mu|}
#       * f_mu(q4, q5)^{2a - N - m - 1}
#       * q4^{||mu^t||^2 / 2} * q5^{-||mu||^2 / 2}
#       * \tilde Z_mu(q4, q5) * \tilde Z_{mu^t}(q5^{-1}, q4^{-1}),
# with f_mu(q4, q5) = (-1)^{|mu|} q4^{||mu^t||^2/2} q5^{||mu||^2/2}.
#
# Writing p = 2a - N - m - 1 and collecting the q-powers, the non-class part is
#     (-1)^{|mu|(p+1)} * r4^{(p+1) ||mu^t||^2} * r5^{(p-1) ||mu||^2}
#         * \tilde Z_mu(q4, q5) * \tilde Z_{mu^t}(q5^{-1}, q4^{-1}).
# The two signs combine the (-1)^{|mu|} of -Q^{[B_a]} and the (-1)^{|mu| p} of f^p.
# The class contribution |mu| [B_a] is handled by the caller.
# ------------------------------------------------------------------------------
function _gp3_one_leg(mu::Vector{Int}, a::Int, N::Int, m::Int, r4, r5, K)
  mut  = _gp3_transpose(mu)
  s    = _gp3_abs(mu)
  nn   = _gp3_sqnorm(mu)       # ||mu||^2
  nnt  = _gp3_sqnorm(mut)      # ||mu^t||^2
  p    = 2 * a - N - m - 1

  sgn   = iseven(s * (p + 1)) ? 1 : -1
  coeff = sgn * _gp3_mono(K, r4, r5, (p + 1) * nnt, (p - 1) * nn)

  coeff *= _gp3_tildeZ(mu,  (2, 0), (0, 2),  r4, r5, K)   # \tilde Z_mu(q4, q5)
  coeff *= _gp3_tildeZ(mut, (0, -2), (-2, 0), r4, r5, K)  # \tilde Z_{mu^t}(q5^{-1}, q4^{-1})
  return coeff
end


# ------------------------------------------------------------------------------
# The character A_{mu,nu} appearing in  R_{mu,nu}(q4,q5;Q) = Exp(Q A_{mu,nu}) * (...),
# evaluated under the k-th Adams operation psi_k (q_i -> q_i^k).
#
#   A_{mu,nu} =  q4 q5^{-1} / ((1-q4)(1-q5^{-1}))
#              + ( sum_{i=1}^{l(nu)} q4^i q5^{-l(mu)-1} (q5^{nu_i} - 1) ) / (1 - q5^{-1})
#              + ( sum_{j=1}^{l(mu)} q5^{-j} q4^{l(nu)+1} (q4^{-mu_j} - 1) ) / (1 - q4)
# ------------------------------------------------------------------------------
function _gp3_R_character(mu::Vector{Int}, nu::Vector{Int}, k::Int, r4, r5, K)
  lmu = length(mu)
  lnu = length(nu)
  oneK = K(1)
  # psi_k(q4^a q5^b) = r4^{2 k a} r5^{2 k b}.
  mon = (a, b) -> _gp3_mono(K, r4, r5, 2 * k * a, 2 * k * b)

  q4   = mon(1, 0)
  q5iv = mon(0, -1)

  res = mon(1, -1) // ((oneK - q4) * (oneK - q5iv))
  for i in 1:lnu
    res += (mon(i, -lmu - 1) * (mon(0, nu[i]) - oneK)) // (oneK - q5iv)
  end
  for j in 1:lmu
    res += (mon(lnu + 1, -j) * (mon(-mu[j], 0) - oneK)) // (oneK - q4)
  end
  return res
end


# ------------------------------------------------------------------------------
# Coefficients [Q^0], [Q^1], ..., [Q^deg] of the series
#
#     R_{mu,nu}(q4, q5; Q)
#       = Exp( Q A_{mu,nu} )
#         * prod_{i=1..l(nu), j=1..l(mu)} (1 - Q q4^i q5^{-j}) / (1 - Q q4^{i-mu_j} q5^{-j+nu_i}).
#
# `Exp` is the plethystic exponential in Q, i.e.
#     Exp(Q A) = exp( sum_{k>=1} (Q^k / k) psi_k(A) ).
# We exponentiate the truncated power sums and multiply by the (rational) finite
# correction expanded in Q.
# ------------------------------------------------------------------------------
function _gp3_R_coeffs(mu::Vector{Int}, nu::Vector{Int}, deg::Int, r4, r5, K)
  deg == 0 && return elem_type(K)[K(1)]

  # log-coefficients L[n] = (1/n) psi_n(A_{mu,nu}),  n = 1..deg.
  L = elem_type(K)[K(0) for _ in 0:deg]
  for n in 1:deg
    L[n + 1] = _gp3_R_character(mu, nu, n, r4, r5, K) // K(n)
  end

  # E = exp(sum_n L[n] Q^n):  E[0] = 1,  E[k] = (1/k) sum_{j=1}^k j L[j] E[k-j].
  E = elem_type(K)[K(0) for _ in 0:deg]
  E[1] = K(1)
  for k in 1:deg
    acc = K(0)
    for j in 1:k
      acc += K(j) * L[j + 1] * E[k - j + 1]
    end
    E[k + 1] = acc // K(k)
  end

  # Finite correction prod (1 - Q alpha_{ij}) / (1 - Q beta_{ij}), expanded in Q.
  fin = elem_type(K)[K(0) for _ in 0:deg]
  fin[1] = K(1)
  for i in 1:length(nu), j in 1:length(mu)
    alpha = _gp3_mono(K, r4, r5, 2 * i, -2 * j)                          # q4^i q5^{-j}
    beta  = _gp3_mono(K, r4, r5, 2 * (i - mu[j]), 2 * (nu[i] - j))        # q4^{i-mu_j} q5^{nu_i-j}
    # Multiply by 1/(1 - beta Q):  tmp[k] = sum_{p=0}^k beta^{k-p} fin[p].
    tmp = elem_type(K)[K(0) for _ in 0:deg]
    for k in 0:deg
      s = K(0)
      for pp in 0:k
        s += beta^(k - pp) * fin[pp + 1]
      end
      tmp[k + 1] = s
    end
    # Multiply by (1 - alpha Q).
    nxt = copy(tmp)
    for k in 1:deg
      nxt[k + 1] -= alpha * tmp[k]
    end
    fin = nxt
  end

  # R = E * fin, truncated to degree deg.
  out = elem_type(K)[K(0) for _ in 0:deg]
  for k in 0:deg
    s = K(0)
    for p in 0:k
      s += E[p + 1] * fin[k - p + 1]
    end
    out[k + 1] = s
  end
  return out
end


# ------------------------------------------------------------------------------
# Bounding box of the enumeration.
#
# The total class of a term is  sum_a |mu_a| [B_a]  +  sum_{a<b} n_{ab} sum_{c=a}^{b-1} [F_c].
# To capture every term with class <= beta we bound, for each generator [B_a] and
# [F_c], its multiplicity over all nonnegative integer representations of beta.
# This is the lattice-point bounding box of  { x >= 0 : M x = beta },  where the
# columns of M are the H_2-classes of the generators.
# ------------------------------------------------------------------------------
function _gp3_preimage_box(B_cols::Vector{NTuple{r,Int}}, F_cols::Vector{NTuple{r,Int}},
                          beta_t::NTuple{r,Int}) where {r}
  N  = length(B_cols)
  nf = length(F_cols)
  Lgen = N + nf

  M = zeros(Int, r, Lgen)
  for i in 1:N, k in 1:r
    M[k, i] = B_cols[i][k]
  end
  for j in 1:nf, k in 1:r
    M[k, N + j] = F_cols[j][k]
  end

  rhs = Int[beta_t[k] for k in 1:r]
  Aineq = zeros(Int, Lgen, Lgen)
  for i in 1:Lgen
    Aineq[i, i] = -1                  # encodes x_i >= 0
  end

  P = polyhedron((Aineq, zeros(Int, Lgen)), (M, rhs))
  @req is_bounded(P) "preimage polytope of beta is unbounded; the gauge series cannot be truncated."
  @req dim(P) >= 0 "beta has no nonnegative representation in the B/F generators."

  # The bounding box maximizes each coordinate over the preimage polytope.  Each
  # maximum is attained at a vertex, so we only inspect the (few) vertices rather
  # than enumerating every lattice point.  Flooring the rational coordinate maxima
  # keeps a valid integer upper bound on every nonnegative monoid preimage of beta.
  verts = vertices(P)
  box = zeros(Int, Lgen)
  for i in 1:Lgen
    box[i] = Int(floor(maximum(v[i] for v in verts)))
  end
  dmax = ntuple(a -> box[a], N)
  emax = nf > 0 ? ntuple(c -> box[N + c], nf) : ()
  return dmax, emax
end


# ------------------------------------------------------------------------------
# Sparse truncated product of two {class -> coeff} dictionaries, dropping any
# monomial whose class exceeds `bound` componentwise.
# ------------------------------------------------------------------------------
function _gp3_dict_mul(a::Dict{NTuple{L,Int},T}, b::Dict{NTuple{L,Int},T},
                      bound::NTuple{L,Int}, K) where {L,T}
  out = Dict{NTuple{L,Int},T}()
  for (ea, ca) in a, (eb, cb) in b
    ev = ntuple(i -> ea[i] + eb[i], L)
    any(ev[i] > bound[i] for i in 1:L) && continue
    out[ev] = get(out, ev, K(0)) + ca * cb
  end
  for k in collect(keys(out))
    iszero(out[k]) && delete!(out, k)
  end
  return out
end


# ------------------------------------------------------------------------------
# Build the truncated partition function Z_{N,m} in *monoid coordinates*.
#
# A term is indexed by the multiplicities of the generators, i.e. by an
# L = 2N-1 tuple  (d_1, ..., d_N, e_1, ..., e_{N-1})  recording how many copies of
# [B_1], ..., [B_N], [F_1], ..., [F_{N-1}] it uses.  These coordinates are always
# nonnegative and bounded by `gamma` (the preimage box), so the componentwise
# truncation is unambiguous -- unlike the H_2 coordinates, in which effective
# classes can have negative entries.  We push the result forward to H_2 only at
# the end (in the public function), which is valid because the pushforward is a
# ring homomorphism and hence commutes with the logarithm.
# ------------------------------------------------------------------------------
function _gp3_build_Z(N::Int, m::Int, gamma::NTuple{L,Int}, r4, r5, K) where {L}
  dmax = ntuple(a -> gamma[a], N)        # bound on |mu_a| = d_a
  part_lists = [_gp3_partitions_upto(dmax[a]) for a in 1:N]

  Z = Dict{NTuple{L,Int}, elem_type(K)}()

  for mus in Iterators.product(part_lists...)
    # One-leg prefactor; its base monoid class is (|mu_1|, ..., |mu_N|, 0, ..., 0).
    leg = K(1)
    for a in 1:N
      leg *= _gp3_one_leg(mus[a], a, N, m, r4, r5, K)
    end
    iszero(leg) && continue
    base = ntuple(i -> i <= N ? _gp3_abs(mus[i]) : 0, L)

    acc = Dict{NTuple{L,Int}, elem_type(K)}()
    acc[base] = leg

    # Multiply in the pair factors  R_{mu_a^t, mu_b}(..; q5 P) R_{mu_a^t, mu_b}(..; q4^{-1} P),
    # P = prod_{c=a}^{b-1} Q^{[F_c]}.  At total fibre power n the coefficient is
    #     sum_{n1+n2=n} R_{n1} q5^{n1} R_{n2} q4^{-n2},
    # and the monoid class adds n to each of the slots e_a, ..., e_{b-1}.
    for a in 1:N-1, b in a+1:N
      md = minimum(gamma[N + c] for c in a:b-1)
      md == 0 && continue

      Rc = _gp3_R_coeffs(_gp3_transpose(mus[a]), mus[b], md, r4, r5, K)

      pser = Dict{NTuple{L,Int}, elem_type(K)}()
      for n in 0:md
        c = K(0)
        for n1 in 0:n
          n2 = n - n1
          c += Rc[n1 + 1] * _gp3_mono(K, r4, r5, 0, 2 * n1) *
               Rc[n2 + 1] * _gp3_mono(K, r4, r5, -2 * n2, 0)
        end
        iszero(c) && continue
        cls = ntuple(i -> (N + a <= i <= N + b - 1) ? n : 0, L)
        pser[cls] = c
      end

      acc = _gp3_dict_mul(acc, pser, gamma, K)
      isempty(acc) && break
    end

    for (cls, c) in acc
      iszero(c) && continue
      Z[cls] = get(Z, cls, K(0)) + c
    end
  end

  for k in collect(keys(Z))
    iszero(Z[k]) && delete!(Z, k)
  end
  return Z
end


# ------------------------------------------------------------------------------
# Connected free energy = ordinary logarithm of the (truncated) series Z.
# With A = Z - 1 (the class-0 coefficient of Z must be 1) we compute
#     log(1 + A) = sum_{n>=1} (-1)^{n-1}/n  A^n,
# truncated so that every class stays <= bound.
# ------------------------------------------------------------------------------
function _gp3_log(Z::Dict{NTuple{L,Int},T}, bound::NTuple{L,Int}, K) where {L,T}
  zero_key = ntuple(_ -> 0, L)
  @assert isone(get(Z, zero_key, K(0))) "Z must have constant term 1."

  A = Dict{NTuple{L,Int},T}()
  for (ev, c) in Z
    ev == zero_key && continue
    A[ev] = c
  end

  out = Dict{NTuple{L,Int},T}()
  isempty(A) && return out

  cur = copy(A)                 # cur holds A^n
  maxpow = sum(bound)
  for n in 1:maxpow
    factor = (iseven(n) ? K(-1) : K(1)) // K(n)   # (-1)^{n-1} / n
    for (ev, c) in cur
      out[ev] = get(out, ev, K(0)) + factor * c
    end
    n == maxpow && break
    cur = _gp3_dict_mul(cur, A, bound, K)
    isempty(cur) && break
  end

  for k in collect(keys(out))
    iszero(out[k]) && delete!(out, k)
  end
  return out
end


# ------------------------------------------------------------------------------
# Substitute  r4 = exp(u t4 / 2),  r5 = exp(u t5 / 2)  into an element of K and
# expand it as a Laurent series in u through degree `max_deg`, returned as an
# element of the fraction field of `parent(u)` (with t4, t5 the user's variables).
#
# We expand over the auxiliary ring Rt = QQ[T4, T5], divide num/den as power
# series in u (after factoring out the u-valuations), and finally map T4 -> t4,
# T5 -> t5 into the user ring.  This mirrors the genus expansion used elsewhere.
# ------------------------------------------------------------------------------
function _gp3_to_user(coeff, t4, t5, u, max_deg::Int)
  Fu = fraction_field(parent(u))
  iszero(coeff) && return zero(Fu)

  num_p = numerator(coeff)
  den_p = denominator(coeff)

  Rt, (T4, T5) = polynomial_ring(QQ, ["T4", "T5"])

  function _maxexp(p)
    a = 0; b = 0
    for ev in AbstractAlgebra.exponent_vectors(p)
      a = max(a, ev[1]); b = max(b, ev[2])
    end
    return (a, b)
  end
  (na, nb) = _maxexp(num_p)
  (da, db) = _maxexp(den_p)

  # Internal u-precision.  Each factor (1 - r_*^a) contributes u-valuation <= 1,
  # so the u-valuation of num/den is bounded by their total r-degree.
  D = max(max_deg, 0) + max(na + nb, da + db, 0) + 6

  _zerov() = elem_type(Rt)[zero(Rt) for _ in 0:D]

  function _mul(a::Vector, b::Vector)
    out = _zerov()
    la = length(a); lb = length(b)
    for i in 1:la
      iszero(a[i]) && continue
      jmax = min(lb, D + 2 - i)
      for j in 1:jmax
        iszero(b[j]) && continue
        out[i + j - 1] += a[i] * b[j]
      end
    end
    return out
  end

  # exp(u * tv / 2) truncated to u-degree D: coefficient of u^k is tv^k / (2^k k!).
  function _exp_half(tv)
    out = _zerov()
    out[1] = one(Rt)
    cp = one(Rt)
    c = QQ(1)
    for k in 1:D
      cp *= tv
      c *= QQ(1, 2 * k)
      out[k + 1] = c * cp
    end
    return out
  end
  r4s = _exp_half(T4)
  r5s = _exp_half(T5)

  function _eval(p, ma::Int, mb::Int)
    p4 = Vector{Vector{elem_type(Rt)}}(undef, ma + 1)
    p4[1] = _zerov(); p4[1][1] = one(Rt)
    for k in 1:ma
      p4[k + 1] = _mul(p4[k], r4s)
    end
    p5 = Vector{Vector{elem_type(Rt)}}(undef, mb + 1)
    p5[1] = _zerov(); p5[1][1] = one(Rt)
    for k in 1:mb
      p5[k + 1] = _mul(p5[k], r5s)
    end
    res = _zerov()
    for (c, ev) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
      term = _mul(p4[ev[1] + 1], p5[ev[2] + 1])
      for i in 1:D + 1
        iszero(term[i]) && continue
        res[i] += c * term[i]
      end
    end
    return res
  end

  num_u = _eval(num_p, na, nb)
  den_u = _eval(den_p, da, db)

  vnum = findfirst(!iszero, num_u)
  vden = findfirst(!iszero, den_u)
  @assert vden !== nothing "denominator vanishes identically after substitution."
  vnum === nothing && return zero(Fu)

  shift = (vnum - 1) - (vden - 1)        # Laurent valuation of the result
  M = max_deg - shift                    # highest internal u-index needed
  M < 0 && return zero(Fu)

  # Precision guard.  The division below needs the substituted numerator and
  # denominator known through internal u-index M (after factoring out their
  # valuations), i.e. D + 1 - vnum >= M and D + 1 - vden >= M.  A connected free
  # energy is a genus expansion starting no lower than u^-2, so shift >= -2 and
  # D = max_deg + max(num/den r-degree) + 6 always suffices; we assert it so that
  # any shortfall is a loud error rather than a silent under-expansion.
  @assert shift >= -2 "connected free energy has a pole worse than u^-2 (shift = $shift)."
  @assert D + 1 - vnum >= M "u-expansion precision too low for the numerator (need larger D)."
  @assert D + 1 - vden >= M "u-expansion precision too low for the denominator (need larger D)."

  Ft = fraction_field(Rt)
  num_f = elem_type(Ft)[Ft(num_u[vnum + i]) for i in 0:min(M, D + 1 - vnum)]
  den_f = elem_type(Ft)[Ft(den_u[vden + i]) for i in 0:min(M, D + 1 - vden)]

  # Power-series division: q[k] = (num[k] - sum_{j>=1} den[j] q[k-j]) / den[0].
  inv0 = inv(den_f[1])
  q = Vector{elem_type(Ft)}(undef, M + 1)
  for k in 0:M
    s = (k + 1) <= length(num_f) ? num_f[k + 1] : Ft(0)
    for j in 1:k
      (j + 1 > length(den_f)) && break
      s -= den_f[j + 1] * q[k - j + 1]
    end
    q[k + 1] = s * inv0
  end

  # Map the Laurent coefficients into the user ring: T4 -> t4, T5 -> t5, u -> u.
  uF = Fu(u)
  out = zero(Fu)
  for k in 0:M
    iszero(q[k + 1]) && continue
    numc = evaluate(numerator(q[k + 1]),   [t4, t5])
    denc = evaluate(denominator(q[k + 1]), [t4, t5])
    val = Fu(numc) // Fu(denc)
    pw = k + shift
    out += pw >= 0 ? val * uF^pw : val // uF^(-pw)
  end
  return out
end


# ------------------------------------------------------------------------------
# Public entry point.
# ------------------------------------------------------------------------------
@doc raw"""
    gkm_5d_gauge_prediction_3(G, t, u, beta, max_genus)

Conjectural $\Omega_\beta$ for the gauge space $\mathcal{Z}_{N,m}=X_{N,m}\times\mathbb{A}^2$
(`G = gkm_5d_gauge(N, m; equiCY=true)`), following the revised specification in
`GAUGE_latex_specs.tex`.

The result is returned in the fraction field of `parent(u)` as a Laurent series in
$u$ through degree `2*max_genus - 2`, with coefficients rational in `t[4]`, `t[5]`
(the affine-plane weights $\epsilon_4,\epsilon_5$).

!!! note
    The GW-gauge correspondence is conjectured for the equivariantly Calabi-Yau
    setup with $|m| < N$ (`GAUGE_latex_specs.tex`).  The series itself is defined
    for any `m`; this function does not restrict `m`, so it can also be evaluated
    on out-of-range examples.  The prediction depends only on `N`, `m` and the
    curve-class lattice of `G`, not on the torus weights, so it is independent of
    the `equiCY` flag; the equivariant-CY substitution is applied downstream by
    [`get_Omega_beta`](@ref) when comparing against actual invariants.
    Following the `cc_mobius` convention, `u` must be the last generator of
    `parent(u)`.
"""
function gkm_5d_gauge_prediction_3(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)
  @req has_attribute(G, :N) && has_attribute(G, :m) "G must be a gauge graph (gkm_5d_gauge)."
  N = get_attribute(G, :N)::Int
  m = get_attribute(G, :m)::Int
  @req max_genus >= 0 "max_genus must be nonnegative."
  @req length(t) >= 5 "t must contain at least the five equivariant parameters t1..t5."

  # Generator classes in H_2:  [B_a] (sections) and [F_c] (fibres).
  B_classes = [curve_class(G, "v$a", "w$a") for a in 1:N]
  F_classes = [curve_class(G, "v$c", "v$(c+1)") for c in 1:N-1]
  @req parent(beta) == parent(B_classes[1]) "beta must be a curve class of G."

  H2 = parent(beta)
  r = rank(H2)
  beta_t = ntuple(i -> Int(beta[i]), r)
  @req any(x -> x != 0, beta_t) "beta must be nonzero."

  # Work over K = Frac(QQ[r4, r5]) with q4 = r4^2, q5 = r5^2.
  R, (r4, r5) = polynomial_ring(QQ, ["r4", "r5"])
  K = fraction_field(R)

  B_cols = NTuple{r,Int}[ntuple(i -> Int(B_classes[a][i]), r) for a in 1:N]
  F_cols = NTuple{r,Int}[ntuple(i -> Int(F_classes[c][i]), r) for c in 1:N-1]

  # Bounding box of the monoid representations of beta -> truncation order gamma.
  dmax, emax = _gp3_preimage_box(B_cols, F_cols, beta_t)
  L = 2 * N - 1
  gamma = ntuple(i -> i <= N ? dmax[i] : emax[i - N], L)

  Zm = _gp3_build_Z(N, m, gamma, r4, r5, K)
  Fm = _gp3_log(Zm, gamma, K)         # connected free energy in monoid coordinates

  # Push forward to H_2:  monoid class (d, e)  ->  sum_a d_a [B_a] + sum_c e_c [F_c].
  Fh2 = Dict{NTuple{r,Int}, elem_type(K)}()
  for (mc, c) in Fm
    img = ntuple(k -> sum((i <= N ? B_cols[i][k] : F_cols[i - N][k]) * mc[i] for i in 1:L), r)
    Fh2[img] = get(Fh2, img, K(0)) + c
  end

  # Convert the connected free energy at beta and all sub-classes beta/k to the
  # user ring, then read off Omega_beta with the package's BPS/Moebius inversion.
  max_deg = 2 * max_genus - 2
  GW = Dict{CC, Any}()
  for gamma_cc in downward_close_ccs([beta])
    g_t = ntuple(i -> Int(gamma_cc[i]), r)
    GW[gamma_cc] = _gp3_to_user(get(Fh2, g_t, K(0)), t[4], t[5], u, max_deg)
  end

  Omega = cc_mobius(GW)
  return Omega[beta]
end
