###
# Truncated algebra for [q_i] and fractions of [q_i].
#
# For all the truncation functionalities below, we assume the input x or ti to have degree one.
###

function t_sq_bracket_frac(num_factors::Vector, denom_factors::Vector, max_deg::Int64)
  @req length(num_factors) > 0 || length(denom_factors) > 0 "need at least a numerator or denominator"

  R = length(num_factors) > 0 ? parent(num_factors[1]) : parent(denom_factors[1])
  S, x = polynomial_ring(fraction_field(R))

  res = one(S) // 1
  max_num_deg = max_deg + length(denom_factors)

  if length(num_factors) > 0
    res *= multi_mullow([t_sq_bracket(x * f, max_num_deg) for f in num_factors], max_num_deg)
  end

  for d in denom_factors
    res *= t_sq_bracket_inv(x * d, max_num_deg)
  end

  d_res = denominator(res)
  println(typeof(res))
  n_res = numerator(res)
  # println(n_res)
  # println(d_res)
  res = truncate(n_res, max_num_deg+1) // d_res
  # println(res)
  return evaluate(res, one(coefficient_ring(R)))
end

function multi_mullow(xs::Vector, max_deg)
  @req length(xs) >= 1 "Need at lest one factor"
  res = truncate(xs[1], max_deg+1)
  for i in 2:length(xs)
    res = mullow(res, xs[i], max_deg + 1)
  end
  return res
end

function t_sq_bracket(ti, max_deg::Int64)
  return 2 * t_sinh(1//2 * ti, max_deg)
end

function t_sq_bracket_inv(ti, max_deg::Int64)
  return 1//2 * t_sinh_inv(1//2 * ti, max_deg)
end

function t_sinh_inv(x, max_deg::Int64)
  return t_sinh_inv_times_x(x, max_deg + 1) // x
end

function x_over_sinh_x_coefficients(max_exponent::Int64)
  # Computes coefficients of x/sinh(x) = Σ c_n * x^{2n} for n = 0, 1, ..., max_exponent÷2
  # Using the formula: c_n = (2 - 4^n) * B_{2n} / (2n)!
  # where B_{2n} are the Bernoulli numbers.
  num_coeffs = div(max_exponent, 2) + 1
  coeffs = Vector{QQFieldElem}(undef, num_coeffs)

  for n in 0:(num_coeffs - 1)
    coeffs[n + 1] = (2 - ZZRingElem(4)^n) * bernoulli(2*n) // factorial(ZZRingElem(2*n))
  end

  return coeffs
end

function t_sinh_inv_times_x(x, max_deg::Int64)
  c = x_over_sinh_x_coefficients(max_deg)
  max_i = div(max_deg, 2)
  x2 = x^2
  powers = [x2^e for e in 0:max_i]
  return sum(powers .* c[1:length(powers)])
end

function t_sinh(x, max_deg::Int64)
  @req max_deg >= 1 "max_deg must be positive"
  x2 = x^2
  x_recent = x
  res  = x
  fctrl = BigInt(1)
  for i in 1:div(max_deg-1, 2)
    fctrl *= (2*i) * (2*i + 1)
    x_recent *= x2
    res += 1//fctrl * x_recent
  end
  return res
end