
# ================================================================================
# Gauge-theory prediction of \Omega_\beta via the plethystic logarithm
# of the partition function \mathcal Z_{N,m} from the preprint (Section 5).
#
# Implements the construction described in `claude_omega_beta_spec (julia specific).md`.
#
# Conventions:
#   - The monoid class is flattened as beta = (d_1,...,d_N, e_1,...,e_{N-1}).
#   - Exact rational expressions in q_4, q_5 are stored via the auxiliary
#     "square-root" variables r4, r5 with q4 = r4^2, q5 = r5^2, because the
#     formulas from the preprint involve half-integer powers of q4, q5.
#   - After the Mobius/Plog sum the result is substituted into a Laurent series
#     ring in u with coefficients in Q(t4, t5) via
#         r4 = exp(u*t4/2),  r5 = exp(u*t5/2).
# ================================================================================

# ---------- Partition utilities ----------

function _part_transpose(mu::Vector{Int})
  isempty(mu) && return Int[]
  mt = zeros(Int, mu[1])
  for x in mu
    for i in 1:x
      mt[i] += 1
    end
  end
  return mt
end

_part_size(mu) = sum(mu; init=0)
_part_norm_sq(mu) = sum(x^2 for x in mu; init=0)

function _partitions_of(n::Int)
  n == 0 && return Vector{Vector{Int}}([Int[]])
  result = Vector{Vector{Int}}()
  _partitions_helper!(result, Int[], n, n)
  return result
end

function _partitions_helper!(result, prefix, remaining, maxpart)
  if remaining == 0
    push!(result, copy(prefix))
    return
  end
  for p in min(remaining, maxpart):-1:1
    push!(prefix, p)
    _partitions_helper!(result, prefix, remaining - p, p)
    pop!(prefix)
  end
end

function _moebius_mu(n::Int)
  n == 1 && return 1
  result = 1
  nn = n
  p = 2
  while p * p <= nn
    if nn % p == 0
      nn = div(nn, p)
      nn % p == 0 && return 0
      result = -result
    end
    p += 1
  end
  nn > 1 && (result = -result)
  return result
end

# ---------- Elementary monomials in the ring K = Frac(Q[r4, r5]) ----------

# Returns r4^a4 * r5^a5 as an element of K.
function _rmono(K, r4, r5, a4::Int, a5::Int)
  return K(r4)^a4 * K(r5)^a5
end

# ---------- One-leg factors and building blocks ----------

# \tilde Z_\mu(t, q) = \prod_{(i,j) \in \mu} (1 - t^(mu_j^t - i + 1) q^(mu_i - j))^{-1}.
#
# `t_exp` and `q_exp` are (a4, a5) exponent pairs such that
#     t = r4^a4 * r5^a5 ,   q = r4^b4 * r5^b5.
function _tilde_Z(mu, t_exp::Tuple{Int,Int}, q_exp::Tuple{Int,Int}, r4, r5, K)
  result = K(1)
  mut = _part_transpose(mu)
  for i in 1:length(mu)
    for j in 1:mu[i]
      tp = mut[j] - i + 1
      qp = mu[i] - j
      a4 = tp * t_exp[1] + qp * q_exp[1]
      a5 = tp * t_exp[2] + qp * q_exp[2]
      factor = K(1) - _rmono(K, r4, r5, a4, a5)
      result = result // factor
    end
  end
  return result
end

# One-leg weight for the site `a` (1-indexed):
#
#     (-Q_a)^{|\mu|} * f_\mu(q4, q5^{-1})^{-(N+m-2a+1)}
#                   * q4^{\|\mu^t\|^2/2} * q5^{-\|\mu\|^2/2}
#                   * \tilde Z_\mu(q4, q5^{-1})
#                   * \tilde Z_{\mu^t}(q5^{-1}, q4).
#
# Returns the above without the (-Q_a)^|\mu| sign (the sign is applied by the caller).
# `use_formula_as_written = true` uses the factor as written in the preprint,
# which may contain a typo: the squared-norms ||mu^t||^2/2 and ||mu||^2/2 are
# suspicious. Setting the flag to false switches to ||mu^t||/2 and ||mu||/2,
# which is the plausible intended version.
function _one_leg(mu, a::Int, N::Int, m::Int, r4, r5, K;
                  use_formula_as_written::Bool)
  mut = _part_transpose(mu)
  nsq  = _part_norm_sq(mu)    # ||mu||  = sum mu_i^2
  nsqt = _part_norm_sq(mut)   # ||mu^t||
  power = -(N + m - 2*a + 1)

  # f_mu(q4, q5^{-1}) = (-1)^{|mu|} * q4^{||mu^t||/2} * (q5^{-1})^{-||mu||/2}
  #                  = (-1)^{|mu|} * r4^{||mu^t||}   * r5^{||mu||}.
  # So f^power contributes sign (-1)^{|mu|*power} and r4, r5 exponents
  #     (nsqt*power, nsq*power).
  sz = _part_size(mu)
  sgn = iseven(sz * power) ? 1 : -1

  if use_formula_as_written
    # q4^{||mu^t||^2/2} * q5^{-||mu||^2/2} = r4^{||mu^t||^2} * r5^{-||mu||^2}.
    extra4 = nsqt^2
    extra5 = -nsq^2
  else
    # q4^{||mu^t||/2} * q5^{-||mu||/2} = r4^{||mu^t||} * r5^{-||mu||}.
    extra4 = nsqt
    extra5 = -nsq
  end

  tot4 = nsqt * power + extra4
  tot5 = nsq  * power + extra5

  coeff = sgn * _rmono(K, r4, r5, tot4, tot5)
  coeff *= _tilde_Z(mu,  (2, 0), (0, -2), r4, r5, K)  # \tilde Z_\mu(q4, q5^{-1})
  coeff *= _tilde_Z(mut, (0, -2), (2, 0), r4, r5, K)  # \tilde Z_{\mu^t}(q5^{-1}, q4)
  return coeff
end

# A_{\lam,\kap}(T, Q) evaluated at (T, Q) = (q4^n, (q5^{-1})^n), i.e. under the
# n-th Adams operation. Returns an element of K.
function _A_expr(lam::Vector{Int}, kap::Vector{Int}, n::Int, r4, r5, K)
  ll = length(lam)
  lk = length(kap)
  # T = q4^n = r4^(2n),  Q = (q5^{-1})^n = r5^(-2n).
  # So T^a * Q^b  =  r4^(2n*a) * r5^(-2n*b).
  e = (a, b) -> _rmono(K, r4, r5, 2*n*a, -2*n*b)
  oneK = K(1)
  T1 = e(1, 0)
  Q1 = e(0, 1)
  result = (T1 * Q1) // ((oneK - T1) * (oneK - Q1))
  for i in 1:ll
    # T^i * Q^(lk+1) * (Q^{-lam_i} - 1) / (1 - Q).
    result += (e(i, lk + 1) * (e(0, -lam[i]) - oneK)) // (oneK - Q1)
  end
  for j in 1:lk
    # Q^j * T^(ll+1) * (T^{-kap_j} - 1) / (1 - T).
    result += (e(0, j) * e(ll + 1, 0) * (e(-kap[j], 0) - oneK)) // (oneK - T1)
  end
  return result
end

# Expand R_{\lam,\kap}(x; q4, q5^{-1}) as a polynomial in x up to degree r_max,
# returning the coefficient list [ [x^0], [x^1], ..., [x^r_max] ] of elements of K.
function _expand_R(lam::Vector{Int}, kap::Vector{Int}, r_max::Int, r4, r5, K)
  r_max == 0 && return elem_type(K)[K(1)]

  # P[k+1] = [x^k] of  \sum_{n=1}^{r_max} (1/n) A(t^n, q^n) * x^n.
  P = elem_type(K)[K(0) for _ in 0:r_max]
  for n in 1:r_max
    P[n+1] = _A_expr(lam, kap, n, r4, r5, K) // K(n)
  end

  # expP = exp(P) truncated to x^r_max.  Use
  #     expP[0] = 1,   expP[k] = (1/k) sum_{j=1}^{k} j * P[j] * expP[k-j].
  expP = elem_type(K)[K(0) for _ in 0:r_max]
  expP[1] = K(1)
  for k in 1:r_max
    s = K(0)
    for j in 1:k
      s += K(j) * P[j+1] * expP[k-j+1]
    end
    expP[k+1] = s // K(k)
  end

  # Finite correction product:
  #   \prod_{i=1}^{ll} \prod_{j=1}^{lk}
  #       (1 - x * alpha_{ij}) / (1 - x * beta_{ij})
  # with alpha_{ij} = t^i q^j          = r4^{ 2 i } r5^{-2 j },
  #      beta_{ij}  = t^{i-kap_j} q^{j-lam_i}
  #                                    = r4^{ 2(i-kap_j)} r5^{-2(j-lam_i)}.
  fin = elem_type(K)[K(0) for _ in 0:r_max]
  fin[1] = K(1)
  ll = length(lam)
  lk = length(kap)
  for i in 1:ll, j in 1:lk
    alpha = _rmono(K, r4, r5,  2*i,            -2*j)
    beta  = _rmono(K, r4, r5,  2*(i - kap[j]), -2*(j - lam[i]))
    # Multiply fin by 1/(1 - beta x): new_fin[k] = sum_{i=0}^{k} beta^{k-i} * fin[i].
    tmp = elem_type(K)[K(0) for _ in 0:r_max]
    for k in 0:r_max
      s = K(0)
      for i2 in 0:k
        s += beta^(k - i2) * fin[i2 + 1]
      end
      tmp[k+1] = s
    end
    # Multiply by (1 - alpha x).
    new_fin = copy(tmp)
    for k in 1:r_max
      new_fin[k+1] -= alpha * tmp[k]
    end
    fin = new_fin
  end

  # Final: expP * fin truncated to x^r_max.
  out = elem_type(K)[K(0) for _ in 0:r_max]
  for k in 0:r_max
    s = K(0)
    for i in 0:k
      s += expP[i+1] * fin[k - i + 1]
    end
    out[k+1] = s
  end
  return out
end

# ---------- Building the truncated Z_{N,m} ----------

# Multiply two sparse dicts representing truncated multivariate series, keeping
# only monomials componentwise <= bound.
function _dict_mul(a::Dict{NTuple{L,Int}, T}, b::Dict{NTuple{L,Int}, T},
                   bound::NTuple{L,Int}, K) where {L, T}
  result = Dict{NTuple{L,Int}, T}()
  for (ea, ca) in a, (eb, cb) in b
    ev = ntuple(i -> ea[i] + eb[i], L)
    any(ev[i] > bound[i] for i in 1:L) && continue
    result[ev] = get(result, ev, K(0)) + ca * cb
  end
  for k in collect(keys(result))
    iszero(result[k]) && delete!(result, k)
  end
  return result
end

# Compute the truncated Z_{N,m}(Q, \tilde Q) as a sparse dictionary
#     exponent tuple (d_1,...,d_N, e_1,...,e_{N-1}) => coefficient in K,
# keeping only monomials componentwise <= gamma.
function _compute_Z_dict(N::Int, m::Int, gamma::NTuple{L,Int},
                         r4, r5, K;
                         use_formula_as_written::Bool) where {L}
  @assert L == 2*N - 1 "beta must have length 2N - 1 = $(2*N-1), got $L"
  d = ntuple(i -> gamma[i],     N)
  e = N > 1 ? ntuple(i -> gamma[N+i], N-1) : ()

  Z = Dict{NTuple{L,Int}, elem_type(K)}()

  part_lists = [reduce(vcat, (_partitions_of(k) for k in 0:d[a])) for a in 1:N]

  for mus in Iterators.product(part_lists...)
    sizes = ntuple(a -> _part_size(mus[a]), N)

    # One-leg prefactor (including the (-Q_a)^{|\mu_a|} sign).
    pref = K(1)
    for a in 1:N
      pref *= _one_leg(mus[a], a, N, m, r4, r5, K;
                       use_formula_as_written = use_formula_as_written)
    end
    if isodd(sum(sizes))
      pref = -pref
    end

    # Precompute pair contributions.  For each ordered pair a < b we get a list
    # of (degree in \tilde Q_{a,b}, coefficient) pairs, obtained from the
    # product of the two R-factors attached to that pair, truncated to degree
    #     r_{a,b} = min(e_a, ..., e_{b-1}).
    pair_contribs = Vector{Vector{Tuple{Int, elem_type(K)}}}()
    if N >= 2
      for a in 1:N-1, b in a+1:N
        r_ab = minimum(e[k] for k in a:b-1)
        if r_ab == 0
          push!(pair_contribs, Tuple{Int, elem_type(K)}[(0, K(1))])
          continue
        end
        lam = _part_transpose(mus[a])
        kap = _part_transpose(mus[b])
        Rcoeffs = _expand_R(lam, kap, r_ab, r4, r5, K)

        # R_{lam,kap}(x1; q4, q5^{-1}) with x1 = \tilde Q_{a,b} * q5
        #   => x1^{n1} = \tilde Q_{a,b}^{n1} * q5^{n1}.
        # R_{lam,kap}(x2; q4, q5^{-1}) with x2 = \tilde Q_{a,b} * q4^{-1}
        #   => x2^{n2} = \tilde Q_{a,b}^{n2} * q4^{-n2}.
        #
        # Combined contribution at total \tilde Q_{a,b}-degree n = n1 + n2 is
        #   sum_{n1+n2=n, 0<=n1,n2<=r_ab}  R[n1+1] * q5^{n1} * R[n2+1] * q4^{-n2}.
        combined = Vector{Tuple{Int, elem_type(K)}}()
        for n in 0:r_ab
          c = K(0)
          for n1 in 0:n
            n2 = n - n1
            n2 > r_ab && continue
            n1 > r_ab && continue
            q5n1     = _rmono(K, r4, r5, 0, 2*n1)
            q4mn2    = _rmono(K, r4, r5, -2*n2, 0)
            c += Rcoeffs[n1+1] * q5n1 * Rcoeffs[n2+1] * q4mn2
          end
          iszero(c) || push!(combined, (n, c))
        end
        push!(pair_contribs, combined)
      end
    end

    # Multiply pair contributions into a sparse dict keyed by the
    # \tilde Q-exponent tuple, pruning monomials beyond the bound.
    EType = NTuple{max(N-1, 0), Int}
    acc = Dict{EType, elem_type(K)}()
    acc[ntuple(_ -> 0, max(N-1, 0))] = K(1)
    pair_idx = 1
    for a in 1:N-1, b in a+1:N
      contribs = pair_contribs[pair_idx]
      pair_idx += 1
      new_acc = Dict{EType, elem_type(K)}()
      for (cur_e, cur_c) in acc
        for (n, c_R) in contribs
          new_e = ntuple(j -> cur_e[j] + ((a <= j <= b-1) ? n : 0), N-1)
          skip = false
          for j in 1:N-1
            if new_e[j] > e[j]
              skip = true
              break
            end
          end
          skip && continue
          nc = cur_c * c_R
          new_acc[new_e] = get(new_acc, new_e, K(0)) + nc
        end
      end
      acc = new_acc
      isempty(acc) && break
    end

    for (ex_e, c_acc) in acc
      iszero(c_acc) && continue
      full_exp = ntuple(i -> i <= N ? sizes[i] : ex_e[i - N], L)
      skip = false
      for i in 1:L
        if full_exp[i] > gamma[i]
          skip = true
          break
        end
      end
      skip && continue
      val = pref * c_acc
      iszero(val) && continue
      Z[full_exp] = get(Z, full_exp, K(0)) + val
    end
  end

  # Clean up zero entries
  for k in collect(keys(Z))
    iszero(Z[k]) && delete!(Z, k)
  end
  return Z
end

# ---------- Plethystic logarithm of the truncated series ----------

function _log_Z(Z_dict::Dict{NTuple{L,Int}, T}, gamma::NTuple{L,Int},
                K) where {L, T}
  zero_exp = ntuple(_ -> 0, L)
  c0 = get(Z_dict, zero_exp, K(0))
  @assert isone(c0) "Z_{N,m} does not have constant term 1 (got $c0)"

  A = Dict{NTuple{L,Int}, T}()
  for (k, v) in Z_dict
    k == zero_exp && continue
    A[k] = v
  end

  result = Dict{NTuple{L,Int}, T}()
  isempty(A) && return result

  A_pow = copy(A)
  max_r = sum(gamma)
  for r in 1:max_r
    sgn = iseven(r) ? K(-1) : K(1)
    scale = sgn // K(r)
    for (ev, c) in A_pow
      result[ev] = get(result, ev, K(0)) + scale * c
    end
    r == max_r && break
    A_pow = _dict_mul(A_pow, A, gamma, K)
    isempty(A_pow) && break
  end

  for k in collect(keys(result))
    iszero(result[k]) && delete!(result, k)
  end
  return result
end

# ---------- Pushforward from the monoid algebra M/~ to H2 ----------

function _pushforward_Z_dict(Z_dict::Dict{NTuple{L,Int}, T},
                             Mmat::Matrix{Int},
                             gamma_h2::NTuple{R,Int},
                             K) where {L, R, T}
  @assert size(Mmat) == (R, L) "expected Mmat to have size ($(R), $(L)), got $(size(Mmat))"

  result = Dict{NTuple{R,Int}, T}()
  for (ev, c) in Z_dict
    img = ntuple(k -> sum(Mmat[k, i] * ev[i] for i in 1:L), R)
    any(img[k] < 0 || img[k] > gamma_h2[k] for k in 1:R) && continue
    result[img] = get(result, img, K(0)) + c
  end

  for k in collect(keys(result))
    iszero(result[k]) && delete!(result, k)
  end
  return result
end

# ---------- Bounding box in M/~ needed to compute a fixed H2-class ----------

function _max_monoid_preimage_bound(Mmat::Matrix{Int}, beta_vec::Vector{Int})
  r, L = size(Mmat)
  @assert length(beta_vec) == r "beta_vec has length $(length(beta_vec)), expected $r"

  A_ineq = zeros(Int, L, L)
  for i in 1:L
    A_ineq[i, i] = -1
  end
  b_ineq = zeros(Int, L)

  P = polyhedron((A_ineq, b_ineq), (Mmat, beta_vec))

  @req is_bounded(P) "preimage polytope under monoid->H2 map is unbounded; ker(pi) meets the nonnegative orthant nontrivially"

  preimages = lattice_points(P)
  @req !isempty(preimages) "no nonnegative integer preimage of beta under the monoid->H2 map"

  gamma = zeros(Int, L)
  for pt in preimages
    for i in 1:L
      gamma[i] = max(gamma[i], Int(pt[i]))
    end
  end

  return ntuple(i -> gamma[i], L)
end

# ---------- Common Möbius extraction step ----------

function _omega_from_log_dict(logZ::Dict{NTuple{L, T1}, T2},
                              beta_t::NTuple{L,Int},
                              r4, r5, K, g_max::Int) where {L, T1, T2}
  nonzero_components = Int[beta_t[i] for i in 1:L if beta_t[i] > 0]
  g = foldl(gcd, nonzero_components)

  total = K(0)
  for k in 1:g
    g % k == 0 || continue
    mu_k = _moebius_mu(k)
    mu_k == 0 && continue
    beta_over_k = ntuple(i -> div(beta_t[i], k), L)
    @assert all(beta_over_k[i] * k == beta_t[i] for i in 1:L)
    coeff_raw = get(logZ, beta_over_k, K(0))
    iszero(coeff_raw) && continue
    coeff_k = _adams_substitute(coeff_raw, k, r4, r5, K)
    total += (K(mu_k) // K(k)) * coeff_k
  end

  return _substitute_exp_and_expand(total, g_max)
end

# ---------- Adams substitution r4 -> r4^k, r5 -> r5^k on a K element ----------

function _adams_substitute(c, k::Int, r4, r5, K)
  k == 1 && return c
  num = numerator(c)
  den = denominator(c)
  num_sub = evaluate(num, [r4^k, r5^k])
  den_sub = evaluate(den, [r4^k, r5^k])
  return K(num_sub) // K(den_sub)
end

# ---------- Substitution r4 = exp(u*t4/2), r5 = exp(u*t5/2) ----------
#
# Strategy: work coefficientwise in u over the polynomial ring Rt = Q[t4, t5]
# (so intermediate series are represented as plain vectors of Rt elements).
# This avoids the expensive gcd/normalization work that happens when the
# coefficient ring is already a fraction field.  After the substitution we
# factor out the u-valuations from numerator and denominator by hand, and only
# then pass to Ft = Frac(Rt) to perform the final series division.

function _substitute_exp_and_expand(coeff, g_max::Int)
  target_prec = 2*g_max - 2  # we want all Laurent terms through u^target_prec

  num_p = numerator(coeff)
  den_p = denominator(coeff)

  Rt, (t4v, t5v) = polynomial_ring(QQ, ["t4", "t5"])

  # Maximum r4, r5 degrees of the numerator / denominator.
  function _max_ev(p)
    a, b = 0, 0
    for ev in AbstractAlgebra.exponent_vectors(p)
      a = max(a, ev[1]); b = max(b, ev[2])
    end
    return (a, b)
  end
  (na, nb) = _max_ev(num_p)
  (da, db) = _max_ev(den_p)

  # Internal u-degree precision.  The u-valuation of each of num(r4,r5) and
  # den(r4,r5) after the substitution is bounded by their total degree in
  # r4, r5 (each factor of the form (1 - r_*^a) contributes u-valuation 1).
  # We need enough room for the shifted quotient up to u^{target_prec}.
  max_r_deg = max(na + nb, da + db, 0)
  D = max(target_prec, 0) + max_r_deg + 5

  # A truncated u-polynomial over Rt is stored as Vector{Rt-elem} of length D+1
  # where index i represents the coefficient of u^{i-1}.
  function _zero_vec()
    return elem_type(Rt)[zero(Rt) for _ in 0:D]
  end

  function _mul_trunc(a::Vector, b::Vector)
    out = _zero_vec()
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

  # exp(u * t / 2) truncated to u-degree D: the coefficient of u^k is
  # t^k / (2^k * k!).
  function _exp_u_half(t_var)
    out = _zero_vec()
    out[1] = one(Rt)
    c_pow = one(Rt)
    coef = QQ(1)
    for k in 1:D
      c_pow *= t_var
      coef *= QQ(1, 2*k)
      out[k+1] = coef * c_pow
    end
    return out
  end

  r4s = _exp_u_half(t4v)
  r5s = _exp_u_half(t5v)

  # Evaluate a QQ[r4, r5] polynomial with r4 -> r4s, r5 -> r5s.
  function _eval_poly(p, max_a::Int, max_b::Int)
    r4p = Vector{Vector{elem_type(Rt)}}(undef, max_a + 1)
    r4p[1] = _zero_vec(); r4p[1][1] = one(Rt)
    for k in 1:max_a
      r4p[k+1] = _mul_trunc(r4p[k], r4s)
    end
    r5p = Vector{Vector{elem_type(Rt)}}(undef, max_b + 1)
    r5p[1] = _zero_vec(); r5p[1][1] = one(Rt)
    for k in 1:max_b
      r5p[k+1] = _mul_trunc(r5p[k], r5s)
    end

    res = _zero_vec()
    for (c, ev) in zip(AbstractAlgebra.coefficients(p),
                       AbstractAlgebra.exponent_vectors(p))
      term = _mul_trunc(r4p[ev[1] + 1], r5p[ev[2] + 1])
      for i in 1:D + 1
        iszero(term[i]) && continue
        res[i] += c * term[i]
      end
    end
    return res
  end

  num_u = _eval_poly(num_p, na, nb)
  den_u = _eval_poly(den_p, da, db)

  # u-valuations (1-indexed positions of the first nonzero entry).
  v_num = findfirst(!iszero, num_u)
  v_den = findfirst(!iszero, den_u)
  @assert v_den !== nothing "denominator vanishes identically after substitution"

  Ft = fraction_field(Rt)

  if v_num === nothing
    Slr, _ = laurent_series_field(Ft, max(target_prec + 1, 1), "u")
    return zero(Slr)
  end

  shift = (v_num - 1) - (v_den - 1)        # final Laurent valuation of the result
  M = target_prec - shift                  # highest u-index of the shifted quotient

  # Precision for the Laurent series ring: each contribution Slr(q[k+1]) * u^(k+shift)
  # inherits precision k+shift+N from the ring's base precision N.  To retain all
  # terms through u^{target_prec} we need 0+shift+N >= target_prec+1, i.e. N >= M+1.
  laurent_prec = max(M + 1, target_prec + 1, 1)
  Slr, u_var = laurent_series_field(Ft, laurent_prec, "u")

  M < 0 && return zero(Slr)

  num_f = elem_type(Ft)[Ft(num_u[v_num + i]) for i in 0:min(M, D + 1 - v_num)]
  den_f = elem_type(Ft)[Ft(den_u[v_den + i]) for i in 0:min(M, D + 1 - v_den)]

  # Power-series division in u: q[k] = (num[k] - sum_{j=1}^{k} den[j] * q[k-j]) / den[0].
  inv_d0 = inv(den_f[1])
  q = Vector{elem_type(Ft)}(undef, M + 1)
  for k in 0:M
    s = (k + 1) <= length(num_f) ? num_f[k + 1] : Ft(0)
    for j in 1:k
      (j + 1 > length(den_f)) && break
      s -= den_f[j + 1] * q[k - j + 1]
    end
    q[k + 1] = s * inv_d0
  end

  result = zero(Slr)
  for k in 0:M
    iszero(q[k + 1]) && continue
    result += Slr(q[k + 1]) * u_var^(k + shift)
  end
  return result
end

# ---------- Public entry point ----------

@doc raw"""
    omega_beta_gauge(N::Int, m::Int, beta, g_max::Int;
                     use_formula_as_written::Bool = true)

Return the Laurent expansion in $u$ of $\Omega_\beta$ for the gauge-theory
partition function $\mathcal Z_{N,m}$ from Section 5 of the preprint, computed
through order $u^{2 g_{\max}-2}$.

The argument `beta` is a tuple or vector of the form
`(d_1, ..., d_N, e_1, ..., e_{N-1})`, interpreted as the monomial
$Q^\beta = \prod_{i=1}^N Q_i^{d_i}\,\prod_{j=1}^{N-1}\tilde Q_j^{e_j}$
in the monoid $M/\sim$ used in the gauge-theory section of the preprint.
It must have nonnegative entries and must be nonzero.

The result is an element of a Laurent series field in $u$ with coefficients in
$\mathbb Q(t_4, t_5)$. Half-integer powers of $q_4, q_5$ that appear in the
intermediate formulas are handled exactly via the auxiliary variables
$r_r = q_r^{1/2}$, and the final substitution is
$r_4 = \exp(u t_4/2),\ r_5 = \exp(u t_5/2)$.

The optional keyword `use_formula_as_written` selects between two versions of
the one-leg weight that differ in the suspicious factor
$q_4^{\|\mu^t\|^2/2} q_5^{-\|\mu\|^2/2}$:

  * `true` (default) -- use the formula **exactly as printed** in the preprint,
    i.e. with the squared norms,
  * `false` -- use the plausible "fixed" version with
    $q_4^{\|\mu^t\|/2} q_5^{-\|\mu\|/2}$ instead.
"""
function omega_beta_gauge(N::Int, m::Int, beta, g_max::Int;
                          use_formula_as_written::Bool = true)

  @req N >= 1 "N must be positive"

  beta_t = Tuple(Int(x) for x in beta)
  L = length(beta_t)
  @req L == 2*N - 1 "length of beta must equal 2N - 1 = $(2*N-1), got $L"
  @req all(x -> x >= 0, beta_t) "beta must have nonnegative entries"
  @req any(x -> x > 0, beta_t) "beta must be nonzero"

  # Work in K = Frac(Q[r4, r5]), with q4 = r4^2, q5 = r5^2.
  R, (r4, r5) = polynomial_ring(QQ, ["r4", "r5"])
  K = fraction_field(R)

  gamma = ntuple(i -> beta_t[i], L)

  Z_dict = _compute_Z_dict(N, m, gamma, r4, r5, K;
                           use_formula_as_written = use_formula_as_written)
  logZ = _log_Z(Z_dict, gamma, K)
  return _omega_from_log_dict(logZ, beta_t, r4, r5, K, g_max)
end

function omega_beta_gauge_h2(N::Int, m::Int, Mmat::Matrix{Int}, beta_h2, g_max::Int;
                             use_formula_as_written::Bool = true)
  @req N >= 1 "N must be positive"

  beta_t = Tuple(Int(x) for x in beta_h2)
  Rdim = length(beta_t)
  @req size(Mmat) == (Rdim, 2*N - 1) "expected Mmat to have size ($(Rdim), $(2*N-1)), got $(size(Mmat))"
  @req all(x -> x >= 0, beta_t) "beta_h2 must have nonnegative entries"
  @req any(x -> x > 0, beta_t) "beta_h2 must be nonzero"

  R, (r4, r5) = polynomial_ring(QQ, ["r4", "r5"])
  K = fraction_field(R)

  beta_vec = [beta_t[i] for i in 1:Rdim]
  gamma_monoid = _max_monoid_preimage_bound(Mmat, beta_vec)

  Z_monoid = _compute_Z_dict(N, m, gamma_monoid, r4, r5, K;
                             use_formula_as_written = use_formula_as_written)
  Z_h2 = _pushforward_Z_dict(Z_monoid, Mmat, beta_t, K)
  logZ_h2 = _log_Z(Z_h2, beta_t, K)

  return _omega_from_log_dict(logZ_h2, beta_t, r4, r5, K, g_max)
end
