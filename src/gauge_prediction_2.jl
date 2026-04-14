# Fresh implementation of the gauge-theory prediction from Section 5 of the
# paper. This file is intentionally independent from `gauge_prediction.jl`.

function _gp2_part_transpose(mu::Vector{Int})
  isempty(mu) && return Int[]
  cols = zeros(Int, mu[1])
  for row in mu
    for j in 1:row
      cols[j] += 1
    end
  end
  return cols
end

_gp2_part_size(mu::Vector{Int}) = sum(mu; init=0)
_gp2_part_quad(mu::Vector{Int}) = sum(x^2 for x in mu; init=0)

function _gp2_partitions(n::Int)
  out = Vector{Vector{Int}}()
  function rec!(prefix::Vector{Int}, remain::Int, mx::Int)
    if remain == 0
      push!(out, copy(prefix))
      return
    end
    for p in min(remain, mx):-1:1
      push!(prefix, p)
      rec!(prefix, remain - p, p)
      pop!(prefix)
    end
  end
  rec!(Int[], n, n)
  return out
end

function _gp2_mobius(n::Int)
  n == 1 && return 1
  x = n
  sign = 1
  p = 2
  while p * p <= x
    if x % p == 0
      x = div(x, p)
      x % p == 0 && return 0
      sign = -sign
    end
    p += 1
  end
  x > 1 && (sign = -sign)
  return sign
end

function _gp2_rpow(K, ρ4, ρ5, a4::Int, a5::Int)
  return K(ρ4)^a4 * K(ρ5)^a5
end

function _gp2_tildeZ(mu::Vector{Int},
                     t_exp::Tuple{Int,Int},
                     q_exp::Tuple{Int,Int},
                     ρ4, ρ5, K)
  ans = K(1)
  muT = _gp2_part_transpose(mu)
  for i in 1:length(mu)
    for j in 1:mu[i]
      ta = muT[j] - i + 1
      qa = mu[i] - j
      e4 = ta * t_exp[1] + qa * q_exp[1]
      e5 = ta * t_exp[2] + qa * q_exp[2]
      ans //= (K(1) - _gp2_rpow(K, ρ4, ρ5, e4, e5))
    end
  end
  return ans
end

function _gp2_leg(mu::Vector{Int}, a::Int, N::Int, m::Int, ρ4, ρ5, K;
                  use_formula_as_written::Bool)
  muT = _gp2_part_transpose(mu)
  sz = _gp2_part_size(mu)
  qmu = _gp2_part_quad(mu)
  qmuT = _gp2_part_quad(muT)
  power = -(N + m - 2*a + 1)

  extra4 = use_formula_as_written ? qmuT^2 : qmuT
  extra5 = use_formula_as_written ? -qmu^2 : -qmu
  sign = iseven(sz * power) ? 1 : -1

  coeff = sign * _gp2_rpow(K, ρ4, ρ5, qmuT * power + extra4, qmu * power + extra5)
  coeff *= _gp2_tildeZ(mu, (2, 0), (0, -2), ρ4, ρ5, K)
  coeff *= _gp2_tildeZ(muT, (0, -2), (2, 0), ρ4, ρ5, K)
  return coeff
end

function _gp2_A_piece(λ::Vector{Int}, ν::Vector{Int}, n::Int, ρ4, ρ5, K)
  ellλ = length(λ)
  ellν = length(ν)
  mon = (a, b) -> _gp2_rpow(K, ρ4, ρ5, 2*n*a, -2*n*b)
  oneK = K(1)
  T = mon(1, 0)
  Q = mon(0, 1)
  out = (T * Q) // ((oneK - T) * (oneK - Q))
  for i in 1:ellλ
    out += (mon(i, ellν + 1) * (mon(0, -λ[i]) - oneK)) // (oneK - Q)
  end
  for j in 1:ellν
    out += (mon(ellλ + 1, 0) * mon(0, j) * (mon(-ν[j], 0) - oneK)) // (oneK - T)
  end
  return out
end

function _gp2_exp_trunc(coeffs::Vector, deg::Int, K)
  out = [K(0) for _ in 0:deg]
  out[1] = K(1)
  for n in 1:deg
    acc = K(0)
    for j in 1:n
      acc += K(j) * coeffs[j+1] * out[n-j+1]
    end
    out[n+1] = acc // K(n)
  end
  return out
end

function _gp2_R_coeffs(λ::Vector{Int}, ν::Vector{Int}, deg::Int, ρ4, ρ5, K)
  base = [K(0) for _ in 0:deg]
  for n in 1:deg
    base[n+1] = _gp2_A_piece(λ, ν, n, ρ4, ρ5, K) // K(n)
  end
  exp_part = _gp2_exp_trunc(base, deg, K)

  corr = [K(0) for _ in 0:deg]
  corr[1] = K(1)
  for i in 1:length(λ), j in 1:length(ν)
    α = _gp2_rpow(K, ρ4, ρ5, 2*i, -2*j)
    β = _gp2_rpow(K, ρ4, ρ5, 2*(i - ν[j]), -2*(j - λ[i]))

    tmp = [K(0) for _ in 0:deg]
    for n in 0:deg
      s = K(0)
      for p in 0:n
        s += corr[p+1] * β^(n-p)
      end
      tmp[n+1] = s
    end

    nxt = copy(tmp)
    for n in 1:deg
      nxt[n+1] -= α * tmp[n]
    end
    corr = nxt
  end

  out = [K(0) for _ in 0:deg]
  for n in 0:deg
    s = K(0)
    for p in 0:n
      s += exp_part[p+1] * corr[n-p+1]
    end
    out[n+1] = s
  end
  return out
end

function _gp2_series_mul(a::Dict{NTuple{L,Int}, T},
                         b::Dict{NTuple{L,Int}, T},
                         bound::NTuple{L,Int},
                         K) where {L, T}
  out = Dict{NTuple{L,Int}, T}()
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

function _gp2_log_series(Z::Dict{NTuple{L,Int}, T}, bound::NTuple{L,Int}, K) where {L, T}
  zero_key = ntuple(_ -> 0, L)
  @assert get(Z, zero_key, K(0)) == K(1)

  A = Dict{NTuple{L,Int}, T}()
  for (ev, c) in Z
    ev == zero_key && continue
    A[ev] = c
  end
  isempty(A) && return Dict{NTuple{L,Int}, T}()

  out = Dict{NTuple{L,Int}, T}()
  cur = copy(A)
  maxpow = sum(bound)
  for n in 1:maxpow
    factor = (iseven(n) ? K(-1) : K(1)) // K(n)
    for (ev, c) in cur
      out[ev] = get(out, ev, K(0)) + factor * c
    end
    n == maxpow && break
    cur = _gp2_series_mul(cur, A, bound, K)
    isempty(cur) && break
  end
  for k in collect(keys(out))
    iszero(out[k]) && delete!(out, k)
  end
  return out
end

function _gp2_adams_on_coeff(c, k::Int, ρ4, ρ5, K)
  k == 1 && return c
  num = evaluate(numerator(c), [ρ4^k, ρ5^k])
  den = evaluate(denominator(c), [ρ4^k, ρ5^k])
  return K(num) // K(den)
end

function _gp2_u_expansion_data(coeff, gmax::Int)
  max_u = 2*gmax - 2
  Rt, (t4v, t5v) = polynomial_ring(QQ, ["t4", "t5"])
  Ft = fraction_field(Rt)

  num = numerator(coeff)
  den = denominator(coeff)

  function max_exponents(poly)
    amax = 0
    bmax = 0
    for ev in AbstractAlgebra.exponent_vectors(poly)
      amax = max(amax, ev[1])
      bmax = max(bmax, ev[2])
    end
    return amax, bmax
  end

  na, nb = max_exponents(num)
  da, db = max_exponents(den)
  cutoff = max(max_u, 0) + max(na + nb, da + db) + 6

  function zero_series()
    return elem_type(Rt)[zero(Rt) for _ in 0:cutoff]
  end

  function mul_trunc(a::Vector, b::Vector)
    out = zero_series()
    for i in eachindex(a)
      iszero(a[i]) && continue
      lim = min(length(b), cutoff + 2 - i)
      for j in 1:lim
        iszero(b[j]) && continue
        out[i + j - 1] += a[i] * b[j]
      end
    end
    return out
  end

  function exp_half(tvar)
    out = zero_series()
    out[1] = one(Rt)
    p = one(Rt)
    c = QQ(1)
    for n in 1:cutoff
      p *= tvar
      c *= QQ(1, 2*n)
      out[n+1] = c * p
    end
    return out
  end

  ρ4ser = exp_half(t4v)
  ρ5ser = exp_half(t5v)

  function eval_poly(poly, amax::Int, bmax::Int)
    p4 = Vector{Vector{elem_type(Rt)}}(undef, amax + 1)
    p5 = Vector{Vector{elem_type(Rt)}}(undef, bmax + 1)
    p4[1] = zero_series(); p4[1][1] = one(Rt)
    p5[1] = zero_series(); p5[1][1] = one(Rt)
    for a in 1:amax
      p4[a+1] = mul_trunc(p4[a], ρ4ser)
    end
    for b in 1:bmax
      p5[b+1] = mul_trunc(p5[b], ρ5ser)
    end
    out = zero_series()
    for (c, ev) in zip(AbstractAlgebra.coefficients(poly), AbstractAlgebra.exponent_vectors(poly))
      term = mul_trunc(p4[ev[1] + 1], p5[ev[2] + 1])
      for i in eachindex(out)
        iszero(term[i]) && continue
        out[i] += c * term[i]
      end
    end
    return out
  end

  numser = eval_poly(num, na, nb)
  denser = eval_poly(den, da, db)
  vnum = findfirst(!iszero, numser)
  vden = findfirst(!iszero, denser)
  @assert vden !== nothing
  if vnum === nothing
    return (Ft, nothing, nothing, elem_type(Ft)[], max(max_u + 1, 1))
  end

  shift = (vnum - 1) - (vden - 1)
  top = max_u - shift
  prec = max(top + 1, max_u + 1, 1)
  top < 0 && return (Ft, shift, top, elem_type(Ft)[], prec)

  ncoefs = [Ft(numser[vnum + i]) for i in 0:min(top, cutoff + 1 - vnum)]
  dcoefs = [Ft(denser[vden + i]) for i in 0:min(top, cutoff + 1 - vden)]
  inv0 = inv(dcoefs[1])
  q = Vector{elem_type(Ft)}(undef, top + 1)
  for n in 0:top
    s = n + 1 <= length(ncoefs) ? ncoefs[n+1] : Ft(0)
    for j in 1:n
      j + 1 > length(dcoefs) && break
      s -= dcoefs[j+1] * q[n-j+1]
    end
    q[n+1] = s * inv0
  end

  return (Ft, shift, top, q, prec)
end

function _gp2_to_u_series(coeff, gmax::Int)
  Ft, shift, top, q, prec = _gp2_u_expansion_data(coeff, gmax)
  L, uL = laurent_series_field(Ft, prec, "u")
  shift === nothing && return zero(L)
  top < 0 && return zero(L)

  out = zero(L)
  for n in 0:top
    iszero(q[n+1]) && continue
    out += L(q[n+1]) * uL^(shift + n)
  end
  return out
end

function _gp2_series_to_user_ring(series, t::Vector, u)
  Fu = fraction_field(parent(u))
  out = zero(Fu)
  uF = Fu(u)
  for n in valuation(series):precision(series)-1
    c = coeff(series, n)
    iszero(c) && continue
    num = evaluate(numerator(c), [t[4], t[5]])
    den = evaluate(denominator(c), [t[4], t[5]])
    cc = Fu(num) // Fu(den)
    out += n >= 0 ? cc * uF^n : cc // uF^(-n)
  end
  return out
end

function _gp2_gauge_generators(G::AbstractGKM_graph, N::Int)
  Qs = [curve_class(G, "v$i", "w$i") for i in 1:N]
  Qts = [curve_class(G, "v$i", "v$(i+1)") for i in 1:N-1]
  return Qs, Qts
end

function _gp2_preimage_bounds(Qs::Vector, Qts::Vector, beta::CC)
  r = rank(parent(beta))
  N = length(Qs)
  L = N + length(Qts)

  cols = Vector{Vector{Int}}(undef, L)
  for i in 1:N
    cols[i] = [Int(Qs[i][k]) for k in 1:r]
  end
  for j in 1:length(Qts)
    cols[N + j] = [Int(Qts[j][k]) for k in 1:r]
  end

  M = zeros(Int, r, L)
  for i in 1:L, k in 1:r
    M[k, i] = cols[i][k]
  end

  rhs = [Int(beta[k]) for k in 1:r]
  A = zeros(Int, L, L)
  for i in 1:L
    A[i, i] = -1
  end
  P = polyhedron((A, zeros(Int, L)), (M, rhs))

  @req is_bounded(P) "nonnegative preimage polytope for beta is unbounded"
  pts = lattice_points(P)
  @req !isempty(pts) "beta has no nonnegative expression in Q/Qt generators"

  γ = zeros(Int, L)
  for p in pts
    for i in 1:L
      γ[i] = max(γ[i], Int(p[i]))
    end
  end

  return ntuple(i -> γ[i], N), ntuple(i -> γ[N + i], length(Qts))
end

function _gp2_add_classes(a::NTuple{L,Int}, b::NTuple{L,Int}) where {L}
  return ntuple(i -> a[i] + b[i], L)
end

function _gp2_scale_class(c::NTuple{L,Int}, n::Int) where {L}
  return ntuple(i -> n * c[i], L)
end

function _gp2_class_tuple(x::CC)
  r = rank(parent(x))
  return ntuple(i -> Int(x[i]), r)
end

function _gp2_build_Z_h2(G::AbstractGKM_graph, beta::CC, N::Int, m::Int, ρ4, ρ5, K;
                         use_formula_as_written::Bool)
  Qs, Qts = _gp2_gauge_generators(G, N)
  beta_t = _gp2_class_tuple(beta)
  r = length(beta_t)
  q_classes = [_gp2_class_tuple(c) for c in Qs]
  qt_classes = [_gp2_class_tuple(c) for c in Qts]
  dmax, emax = _gp2_preimage_bounds(Qs, Qts, beta)

  part_lists = [vcat([Int[]], reduce(vcat, (_gp2_partitions(n) for n in 1:dmax[a]); init=Vector{Vector{Int}}())) for a in 1:N]

  Z = Dict{NTuple{r,Int}, elem_type(K)}()

  for μs in Iterators.product(part_lists...)
    leg_coeff = K(1)
    total_class = ntuple(_ -> 0, r)
    total_boxes = 0
    for a in 1:N
      μ = μs[a]
      s = _gp2_part_size(μ)
      total_boxes += s
      leg_coeff *= _gp2_leg(μ, a, N, m, ρ4, ρ5, K;
                            use_formula_as_written = use_formula_as_written)
      total_class = _gp2_add_classes(total_class, _gp2_scale_class(q_classes[a], s))
    end
    isodd(total_boxes) && (leg_coeff = -leg_coeff)
    any(total_class[i] > beta_t[i] for i in 1:r) && continue

    acc = Dict{NTuple{r,Int}, elem_type(K)}()
    acc[total_class] = leg_coeff

    for a in 1:N-1, b in a+1:N
      maxdeg = minimum(emax[k] for k in a:b-1)
      maxdeg == 0 && continue

      λ = _gp2_part_transpose(μs[a])
      ν = _gp2_part_transpose(μs[b])
      coeffs = _gp2_R_coeffs(λ, ν, maxdeg, ρ4, ρ5, K)
      seg_class = ntuple(i -> begin
        s = 0
        for k in a:b-1
          s += qt_classes[k][i]
        end
        s
      end, r)

      pair_series = Dict{NTuple{r,Int}, elem_type(K)}()
      pair_series[ntuple(_ -> 0, r)] = K(1)
      for n in 1:maxdeg
        c = K(0)
        for n1 in 0:n
          n2 = n - n1
          n1 > maxdeg && continue
          n2 > maxdeg && continue
          c += coeffs[n1+1] * _gp2_rpow(K, ρ4, ρ5, 0, 2*n1) *
               coeffs[n2+1] * _gp2_rpow(K, ρ4, ρ5, -2*n2, 0)
        end
        iszero(c) && continue
        cls = _gp2_scale_class(seg_class, n)
        any(cls[i] < 0 || cls[i] > beta_t[i] for i in 1:r) && continue
        pair_series[cls] = c
      end

      acc = _gp2_series_mul(acc, pair_series, beta_t, K)
      isempty(acc) && break
    end

    for (cls, c) in acc
      Z[cls] = get(Z, cls, K(0)) + c
    end
  end

  for k in collect(keys(Z))
    iszero(Z[k]) && delete!(Z, k)
  end
  return Z
end

function omega_beta_gauge_2(G::AbstractGKM_graph, beta::CC, gmax::Int;
                            use_formula_as_written::Bool = true)
  @req has_attribute(G, :N) && has_attribute(G, :m) "G must be a gauge graph"
  N = get_attribute(G, :N)::Int
  m = get_attribute(G, :m)::Int
  beta_t = _gp2_class_tuple(beta)
  @req any(x -> x > 0, beta_t) "beta must be nonzero"
  @req all(x -> x >= 0, beta_t) "beta must have nonnegative coordinates in the chosen H2 basis"

  S, (ρ4, ρ5) = polynomial_ring(QQ, ["rho4", "rho5"])
  K = fraction_field(S)

  Z = _gp2_build_Z_h2(G, beta, N, m, ρ4, ρ5, K;
                      use_formula_as_written = use_formula_as_written)
  logZ = _gp2_log_series(Z, beta_t, K)

  nz = Int[x for x in beta_t if x > 0]
  g = foldl(gcd, nz)
  total = K(0)
  for k in 1:g
    g % k == 0 || continue
    μk = _gp2_mobius(k)
    μk == 0 && continue
    βk = ntuple(i -> div(beta_t[i], k), length(beta_t))
    @assert all(βk[i] * k == beta_t[i] for i in 1:length(beta_t))
    coeff = get(logZ, βk, K(0))
    iszero(coeff) && continue
    total += (K(μk) // K(k)) * _gp2_adams_on_coeff(coeff, k, ρ4, ρ5, K)
  end

  return _gp2_to_u_series(total, gmax)
end

function gkm_5d_gauge_prediction_2(G::AbstractGKM_graph, t::Vector, u, beta::CC, max_genus::Int64)
  series = omega_beta_gauge_2(G, beta, max_genus; use_formula_as_written=false)
  return _gp2_series_to_user_ring(series, t, u)
end
