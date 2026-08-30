
# Helper: truncate a multivariate polynomial to total degree <= deg
function _truncate(f::MPolyRingElem, deg::Int)
  R = parent(f)
  ctx = MPolyBuildCtx(R)
  for (c, ev) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    if sum(ev) <= deg
      push_term!(ctx, c, ev)
    end
  end
  return finish(ctx)
end

# # Helper: compute exp(x) as a truncated power series in a multivariate polynomial ring
# # (Superseded by the multinomial coefficient approach in _substitute_exp_in_poly)
# function _exp_trunc(x::MPolyRingElem, deg::Int)
#   R = parent(x)
#   result = one(R)
#   x_power = one(R)
#   for k in 1:deg
#     x_power = _truncate(x_power * x, deg)
#     iszero(x_power) && break
#     result += inv(QQ(factorial(ZZ(k)))) * x_power
#   end
#   return result
# end

# Sparse polynomial represented as Dict{Vector{Int}, QQFieldElem}.
# Used during parsing to avoid O(n²) cost of sequential Oscar polynomial additions.
const _TermDict = Dict{Vector{Int}, QQFieldElem}

# a := a + b
function _dict_add!(a::_TermDict, b::_TermDict)
  for (ev, c) in b
    a[ev] = get(a, ev, QQ(0)) + c
  end
  return a
end

# a := a - b
function _dict_sub!(a::_TermDict, b::_TermDict)
  for (ev, c) in b
    a[ev] = get(a, ev, QQ(0)) - c
  end
  return a
end

# Return -a, not changing a.
function _dict_neg(a::_TermDict)
  return _TermDict(ev => -c for (ev, c) in a)
end

# return a * b, not changing a or b
function _dict_mul(a::_TermDict, b::_TermDict)
  result = _TermDict()
  sizehint!(result, length(a) * length(b))
  for (ev1, c1) in a, (ev2, c2) in b
    ev = ev1 .+ ev2
    result[ev] = get(result, ev, QQ(0)) + c1 * c2
  end
  return result
end

function _dict_pow(a::_TermDict, n::Int)
  n == 0 && return _TermDict(zeros(Int, length(first(keys(a)))) => QQ(1))
  # Repeated squaring
  temp = a
  result = nothing
  bits = n
  while bits > 0
    if bits & 1 == 1
      result = result === nothing ? temp : _dict_mul(result, temp)
    end
    bits >>= 1
    bits > 0 && (temp = _dict_mul(temp, temp))
  end
  return result
end

function _dict_to_poly(d::_TermDict, R::MPolyRing)
  ctx = MPolyBuildCtx(R)
  for (ev, c) in d
    iszero(c) && continue
    push_term!(ctx, c, ev)
  end
  return finish(ctx)
end

# Evaluate a Julia expression tree into a _TermDict (iterative, avoids stack overflow)
function _eval_to_dict(root_expr, var_dicts::Dict{Symbol, _TermDict}, n_vars::Int)
  work = Any[(:eval, root_expr)]
  values_stack = _TermDict[]

  while !isempty(work)
    tag, payload = pop!(work)

    if tag === :apply
      if payload isa Tuple{Symbol, Int, Int}
        # (^, 1, exponent)
        _, _, exp_val = payload
        base = pop!(values_stack)
        push!(values_stack, _dict_pow(base, exp_val))
      else
        op, nargs = payload
        args = _TermDict[pop!(values_stack) for _ in 1:nargs]
        reverse!(args)
        if op === :+
          result = args[1]
          for i in 2:nargs
            _dict_add!(result, args[i])
          end
          push!(values_stack, result)
        elseif op === :*
          result = args[1]
          for i in 2:nargs
            result = _dict_mul(result, args[i])
          end
          push!(values_stack, result)
        elseif op === :- && nargs == 1
          push!(values_stack, _dict_neg(args[1]))
        elseif op === :- && nargs == 2
          _dict_sub!(args[1], args[2])
          push!(values_stack, args[1])
        else
          error("Unknown operator $op with $nargs arguments")
        end
      end

    else  # tag === :eval
      expr = payload
      if expr isa Symbol
        push!(values_stack, copy(var_dicts[expr]))
      elseif expr isa Integer
        push!(values_stack, _TermDict(zeros(Int, n_vars) => QQ(expr)))
      elseif expr isa Expr && expr.head == :call
        op = expr.args[1]::Symbol
        operands = expr.args[2:end]
        if op === :^ && length(operands) == 2 && operands[2] isa Integer
          push!(work, (:apply, (:^, 1, operands[2])))
          push!(work, (:eval, operands[1]))
        else
          n = length(operands)
          push!(work, (:apply, (op, n)))
          for i in n:-1:1
            push!(work, (:eval, operands[i]))
          end
        end
      else
        error("Cannot parse expression: $expr (type: $(typeof(expr)))")
      end
    end
  end

  @assert length(values_stack) == 1
  return values_stack[1]
end

# Parse a polynomial expression string into an Oscar ring element
function _parse_poly_expr(s::String, R::MPolyRing, var_dicts::Dict{Symbol, _TermDict}, n_vars::Int)
  expr = Meta.parse(s)
  d = _eval_to_dict(expr, var_dicts, n_vars)
  return _dict_to_poly(d, R)
end

# Helper: compute the power series quotient num/den truncated at total degree deg
function _series_div(num::MPolyRingElem, den::MPolyRingElem, deg::Int)
  R = parent(num)
  c0 = constant_coefficient(den)
  @assert !iszero(c0) "Denominator has zero constant term after exp substitution"

  # den = c0 * (1 + g) where g has no constant term
  inv_c0 = inv(c0)
  g = _truncate(inv_c0 * den - one(R), deg)

  # 1/(1+g) = sum_{k=0}^{deg} (-g)^k, truncated
  inv_1_plus_g = one(R)
  neg_g_power = one(R)
  for k in 1:deg
    neg_g_power = _truncate(neg_g_power * (-g), deg)
    iszero(neg_g_power) && break
    inv_1_plus_g += neg_g_power
  end

  return _truncate(inv_c0 * num * inv_1_plus_g, deg)
end

# Enumerate all exponent vectors in n_vars variables with total degree <= deg
function _enumerate_exps(n_vars::Int, deg::Int)
  result = Vector{Vector{Int}}()
  _enum_exps_helper!(result, zeros(Int, n_vars), 1, deg, n_vars)
  return result
end

function _enum_exps_helper!(result, current, pos, remaining, n_vars)
  if pos > n_vars
    push!(result, copy(current))
    return
  end
  for k in 0:remaining
    current[pos] = k
    _enum_exps_helper!(result, current, pos + 1, remaining - k, n_vars)
  end
  current[pos] = 0
end

# Substitute qi = exp(ti) using direct multinomial coefficient computation.
# The input polynomial must have integer coefficients, despite being over QQ.
# For each monomial c*q^a, exp(a·t) = sum_{|k|<=deg} prod(a_i^k_i/k_i!) * t^k.
# Accumulates contributions using Int128 arithmetic, avoiding Oscar polynomial ops.
function _substitute_exp_in_poly(f::QQMPolyRingElem, t::Vector{QQMPolyRingElem}, deg::Int)
  R = parent(t[1])
  iszero(f) && return zero(R)
  n_vars = length(t)

  # Certify that every fixed-width intermediate below is representable.
  coeff_abs_sum = ZZ(0)
  max_a = 0
  for (c, ev) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    isone(denominator(c)) || error("exponential substitution requires integer coefficients")
    coeff_abs_sum += abs(numerator(c))
    max_a = max(max_a, maximum(ev))
  end
  factorial(ZZ(deg)) <= typemax(Int) || error("factorials do not fit in Int")
  coeff_abs_sum * ZZ(max(1, max_a))^deg <= typemax(Int128) || error("exponential substitution may overflow Int128")

  # Enumerate target exponent vectors
  target_exps = _enumerate_exps(n_vars, deg)
  n_targets = length(target_exps)

  # Precompute factorial products for each target exponent: prod(k_i!)
  fact_denoms = Vector{Int128}(undef, n_targets)
  for (j, k) in enumerate(target_exps)
    fact_denoms[j] = Int128(1)
    for ki in k
      fact_denoms[j] *= Int128(factorial(ki))
    end
  end

  # Precompute power table: pow_table[a+1, k+1] = a^k as Int128
  pow_table = Matrix{Int128}(undef, max_a + 1, deg + 1)
  for a in 0:max_a
    pow_table[a + 1, 1] = Int128(1)
    for k in 1:deg
      pow_table[a + 1, k + 1] = pow_table[a + 1, k] * Int128(a)
    end
  end

  # Accumulate: for each target exp k, sum over monomials of c * prod(a_i^k_i)
  # Use Int128 to avoid GC-heavy ZZ/QQ allocations in the inner loop
  accum = zeros(Int128, n_targets)

  for (c, ev) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    c_int = Int128(numerator(c))
    for j in 1:n_targets
      k = target_exps[j]
      prod_powers = Int128(1)
      for i in 1:n_vars
        prod_powers *= pow_table[ev[i] + 1, k[i] + 1]
      end
      accum[j] += c_int * prod_powers
    end
  end

  # Build Oscar polynomial: coeff = accum[j] / fact_denoms[j]
  ctx = MPolyBuildCtx(R)
  for j in 1:n_targets
    iszero(accum[j]) && continue
    push_term!(ctx, QQ(ZZ(accum[j]), ZZ(fact_denoms[j])), target_exps[j])
  end
  return finish(ctx)
end

# Build the shared variable dictionaries used for parsing
function _make_var_dicts(n_vars::Int)
  function _single_var_dict(idx)
    ev = zeros(Int, n_vars)
    ev[idx] = 1
    return _TermDict(ev => QQ(1))
  end

  e1_dict = _TermDict()
  for i in 1:4; ev = zeros(Int, n_vars); ev[i] = 1; e1_dict[ev] = QQ(1); end

  e2_pairs = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
  e2_dict = _TermDict()
  for (i,j) in e2_pairs; ev = zeros(Int, n_vars); ev[i] = 1; ev[j] = 1; e2_dict[ev] = QQ(1); end

  e3_triples = [(1,2,3),(1,2,4),(1,3,4),(2,3,4)]
  e3_dict = _TermDict()
  for (i,j,k) in e3_triples; ev = zeros(Int, n_vars); ev[i] = 1; ev[j] = 1; ev[k] = 1; e3_dict[ev] = QQ(1); end

  e4_dict = _TermDict([1,1,1,1,0] => QQ(1))

  return Dict{Symbol, _TermDict}(
    :q0 => _single_var_dict(1), :q1 => _single_var_dict(2),
    :q2 => _single_var_dict(3), :q3 => _single_var_dict(4), :q4 => _single_var_dict(5),
    :e1 => e1_dict, :e2 => e2_dict, :e3 => e3_dict, :e4 => e4_dict,
  )
end

@doc raw"""
    load_omega_data(prefix::String, degrees; data_dir::String)

Load Omega data from files `{prefix}_num.dat` and `{prefix}_den.dat`
in `data_dir` and return a `Dict{Int, Q}` mapping each requested degree $d$
to the rational function $Q_d$ in the fraction field of $\mathbb{Q}[q_0,\ldots,q_4]$.

`degrees` can be any collection of positive integers, or a single `Int` $d_{\max}$
(equivalent to `1:d_{\max}`). Degrees exceeding the number of lines in the files are ignored.

The files may use elementary symmetric polynomials `e1,e2,e3,e4` in `q0,...,q3`,
which are automatically substituted.
"""
function load_omega_data(prefix::String, degrees;
                         data_dir::String=joinpath(@__DIR__, "..", "test", "omega_data"))
  num_file = joinpath(data_dir, "$(prefix)_num.dat")
  den_file = joinpath(data_dir, "$(prefix)_den.dat")

  num_lines = readlines(num_file)
  den_lines = readlines(den_file)
  @assert length(num_lines) == length(den_lines) "num and den files must have the same number of lines"

  N = length(num_lines)
  ds = sort(unique([d for d in degrees if 1 <= d <= N]))

  R, _ = polynomial_ring(QQ, ["q0", "q1", "q2", "q3", "q4"])
  n_vars = 5
  var_dicts = _make_var_dicts(n_vars)
  F = fraction_field(R)

  result = Dict{Int, elem_type(F)}()
  for d in ds
    num_poly = _parse_poly_expr(num_lines[d], R, var_dicts, n_vars)
    den_poly = _parse_poly_expr(den_lines[d], R, var_dicts, n_vars)
    result[d] = F(num_poly) // F(den_poly)
  end

  return result
end

# Convenience method: load_omega_data(prefix, d_max) loads degrees 1:d_max
load_omega_data(prefix::String, d_max::Int; kwargs...) =
  load_omega_data(prefix, 1:d_max; kwargs...)

@doc raw"""
    omega_exp_substitution(Qs::Dict{Int}, deg::Int)

Substitute $q_i = \exp(t_i)$ in each rational function in `Qs` and expand
as a power series in $t_0,\ldots,t_4$ truncated at total degree `deg`.
The numerator and denominator polynomials must have integer coefficients.

`Qs` is a `Dict{Int, Q}` as returned by [`load_omega_data`](@ref).

Returns a `Dict{Int, poly}` with the same keys, where each value is a polynomial
in $\mathbb{Q}[t_0,\ldots,t_4]$.
"""
function omega_exp_substitution(Qs::Dict{Int}, deg::Int)
  R, t = polynomial_ring(QQ, ["t0", "t1", "t2", "t3", "t4"])

  result = Dict{Int, QQMPolyRingElem}()
  for (d, Q) in Qs
    num_sub = _substitute_exp_in_poly(numerator(Q), t, deg)
    den_sub = _substitute_exp_in_poly(denominator(Q), t, deg)
    result[d] = _series_div(num_sub, den_sub, deg)
  end

  return result, t
end

@doc raw"""
    omega_u_scaling(polys::Dict{Int}, T::Vector, u)

Substitute $t_i = u \cdot T_i$ in each polynomial.

`polys` is a `Dict{Int, poly}` as returned by [`omega_exp_substitution`](@ref),
`T` is a vector of variables $[T_0,\ldots,T_4]$, and `u` is a variable,
all living in the same polynomial ring.

Returns a `Dict{Int, poly}` with the same keys.
"""
function omega_u_scaling(polys::Dict{Int}, T::Vector, u)
  subst = [u * T[i] for i in 1:length(T)]
  return Dict{Int, eltype(values(polys))}(d => evaluate(p, subst) for (d, p) in polys)
end

@doc raw"""
    omega_series(prefix::String, degrees, T::Vector, u, deg::Int; data_dir::String)

Combined function: load Omega data for the requested degrees, substitute
$q_i = \exp(u \cdot T_i)$, and expand to total degree `deg` in the $T_i$ and $u$.

`degrees` can be any collection of positive integers, or a single `Int` $d_{\max}$
(equivalent to `1:d_{\max}`).

Returns a `Dict{Int, poly}` mapping each degree to the resulting polynomial in the
ring of `T` and `u`.
"""
function omega_series(prefix::String, degrees, T::Vector, u, deg::Int;
                      data_dir::String=joinpath(@__DIR__, "..", "test", "omega_data"))
  Qs = load_omega_data(prefix, degrees; data_dir=data_dir)
  polys, _ = omega_exp_substitution(Qs, deg)
  return omega_u_scaling(polys, T, u)
end
