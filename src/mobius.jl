@doc raw"""
    cc_mobius(GW::Dict{CC, Any}) -> Dict{CC, Any}

# Input:
* `GW::Dict{CC, Any}`: For each curve class $\beta$, this dictionary contains the coefficient of
   $Q^\beta$ in $GW$.

# Output:
* For each curve class `\beta` of the original input, the resulting dictionary contains `\Omega_{\beta}` as defined below.

# Relation between coefficients of $GW$ and $\Omega_\beta$:
For each key $\beta$ of the input, the $Q^\beta$-coefficient of $GW$ is an expression in some variables $\varepsilon_i$ and
the variable $u$. By convention, $u$ is the last generator.
Then $\Omega_\beta$ is an expression in the same variables.

The key relation is:

$$GW = \sum_{\beta\neq = 0}\sum_{k>0} Q^{k\beta}/k \Omega_\beta(\varepsilon_i, ku).$$

!!! note
    If $\beta$ is a key of `GW` and $k>1$ divides $\beta$ then the coefficient of `Q^{\beta/k}$ in $GW$ is assumed to be zero if not in `GW`.

!!! warning
    The zero curve class should not be set in `GW`.
"""
function cc_mobius(GW::Dict{CC, Any})::Dict{CC, Any}

  @req all(beta -> !iszero(beta), keys(GW)) "GW should not contain the zero class."

  res = Dict{CC, Any}()
  # create a copy of input with zeros filled in where necessary.
  GW = cc_mobius_add_zeros(GW)

  # Sort curve classes by GCD (ascending order)
  # This ensures we compute Ω_{β/k} before Ω_β for any k > 1
  sorted_betas = sort(collect(keys(GW)), by=gcd_curve_class)

  # For each curve class β in ascending order of GCD
  for beta in sorted_betas
    # Skip the zero class if it exists
    if iszero(beta)
      continue
    end

    # Solve for Ω_β from the recursion:
    # GW_β = (1/1) * Ω_β(ε_i, u) + Σ_{k>1, k|gcd(β)} (1/k) * Ω_{β/k}(ε_i, k*u)
    # Therefore: Ω_β = GW_β - Σ_{k>1, k|gcd(β)} (1/k) * Ω_{β/k}(ε_i, k*u)
    res[beta] = compute_omega_recursive(beta, GW, res)
  end

  return res
end

"""
    gcd_curve_class(beta::CC) -> Int

Compute the GCD of all entries in the curve class beta.
"""
function gcd_curve_class(beta::CC)
  return gcd([beta[i] for i in 1:ngens(parent(beta))])
end

"""
    divisors(n::ZZRingElem) -> Vector{Int}

Return all positive divisors of the absolute value of n in ascending order.
"""
function divisors(n::ZZRingElem)
  # Work with the absolute value
  n_abs = abs(n)

  if n_abs == 0
    return Int[]
  end

  divs = Int[]
  # Convert to Int for iteration
  n_int = Int(n_abs)

  for i in 1:isqrt(n_int)
    if n_int % i == 0
      push!(divs, i)
      if i != n_int ÷ i
        push!(divs, n_int ÷ i)
      end
    end
  end

  sort!(divs)
  return divs
end

"""
    compute_omega_recursive(beta::CC, GW::Dict{CC, Any}, Omega::Dict{CC, Any}) -> Any

Compute Ω_β by solving the recursion directly:
GW_β = (1/1) * Ω_β(ε_i, u) + Σ_{k>1, k|gcd(β)} (1/k) * Ω_{β/k}(ε_i, k*u)

Rearranging: Ω_β = GW_β - Σ_{k>1, k|gcd(β)} (1/k) * Ω_{β/k}(ε_i, k*u)

The Omega dictionary should already contain Ω_{β/k} for all k > 1 that divide gcd(β).
"""
function compute_omega_recursive(beta::CC, GW::Dict{CC, Any}, Omega::Dict{CC, Any})
  # Start with GW_β
  result = GW[beta]

  # Get the GCD of all entries in beta
  g = gcd_curve_class(beta)

  # Find all divisors of g greater than 1
  divs = divisors(g)

  # Subtract contributions from Ω_{β/k} for k > 1
  for k in divs
    if k == 1
      continue
    end

    # Compute β/k
    beta_div_k = (parent(beta))([beta[i] ÷ k for i in 1:ngens(parent(beta))])

    # Check if Ω_{β/k} has been computed
    if !haskey(Omega, beta_div_k)
      # error("Omega lacks key $beta_div_k")
      continue
    end

    # Subtract (1/k) * Ω_{β/k}(ε_i, k*u)
    result -= (1 // k) * substitute_u(Omega[beta_div_k], k)
  end

  return result
end

"""
    substitute_u(expr, k::Integer) -> Any

Substitute u with k*u in the expression.
The variable u is assumed to be the last generator.
"""
function substitute_u(expr, k::Integer)
  # Get the parent ring of the expression
  R = parent(expr)

  # Get the number of generators
  n = ngens(R)

  if n == 0
    return expr
  end

  # Create the substitution: replace the last generator (u) with k times itself
  # All other generators stay the same
  new_gens = [gen(R, i) for i in 1:(n-1)]
  push!(new_gens, k * gen(R, n))

  # Perform the substitution using evaluate
  return evaluate(expr, new_gens)
end

@doc raw"""
    cc_mobius_invert(Omega::Dict{CC, Any}; max_k::Int=10) -> Dict{CC, Any}

Test function that computes GW from $\Omega$ using the formula:
$$GW = \sum_{\beta \neq 0} \sum_{k>0} \frac{Q^{k\beta}}{k} \Omega_\beta(\varepsilon_i, ku)$$

This is the inverse operation of `cc_mobius` and can be used to verify correctness.

# Input:
* `Omega::Dict{CC, Any}`: Dictionary containing $\Omega_\beta$ for each curve class $\beta$
* `max_k::Int`: Maximum value of k to use in the sum (default: 10)

# Output:
* `Dict{CC, Any}`: Dictionary containing the computed GW invariants
"""
function cc_mobius_invert(Omega::Dict{CC, Any}; max_k::Int=10)
  GW = Dict{CC, Any}()

  # For each curve class β in Omega, compute contributions to GW
  for (beta, omega_beta) in Omega
    # Skip the zero class if it exists
    if iszero(beta)
      continue
    end

    # For each multiplier k, add the contribution to GW_{k*β}
    for k in 1:max_k
      k_beta = k * beta

      # Compute the contribution: (1/k) * Ω_β(ε_i, k*u)
      contribution = (1 // k) * substitute_u(omega_beta, k)

      # Add to GW_{k*β}
      if haskey(GW, k_beta)
        GW[k_beta] += contribution
      else
        GW[k_beta] = contribution
      end
    end
  end

  return GW
end

"""
    cc_mobius_add_zeros(GW::Dict{CC, Any}) -> Dict{CC, Any}

Create a copy of GW and add zero entries for all curve classes β/k where β is in GW and k > 1 divides gcd(β).

This ensures that the divisibility requirement in the docstring of `cc_mobius` is satisfied:
"if β is a key of GW and k>1 divides β, then β/k is also a key of GW."
"""
function cc_mobius_add_zeros(GW::Dict{CC, Any})
  # Start with a copy of the original dictionary
  GW_extended = copy(GW)

  # Get a sample value to determine the zero element
  sample_beta = first(keys(GW))
  sample_value = GW[sample_beta]
  zero_value = zero(sample_value)

  # Keep track of curve classes we need to add
  to_add = Set{CC}()

  # For each curve class in GW, find all its divisors
  for beta in keys(GW)
    if iszero(beta)
      continue
    end

    g = gcd_curve_class(beta)
    divs = divisors(g)

    # For each divisor k > 1, add β/k if it's not already present
    for k in divs
      if k == 1
        continue
      end

      beta_div_k = (parent(beta))([beta[i] ÷ k for i in 1:ngens(parent(beta))])

      if !haskey(GW_extended, beta_div_k) && !(beta_div_k in to_add)
        push!(to_add, beta_div_k)
      end
    end
  end

  # Add all the missing curve classes with zero value
  for beta in to_add
    GW_extended[beta] = zero_value
  end

  return GW_extended
end

# For debugging only.
function cc_mobius_test(GW::Dict{CC, Any})
  Omega = cc_mobius(GW)
  GW = cc_mobius_add_zeros(GW)
  GW_test = cc_mobius_invert(Omega)
  for (beta, GW_beta) in GW
    iszero(beta) && continue
    if GW_test[beta] != GW_beta
      error("Disagreement for $beta.")
    end
  end
  println("Classes agree.")
  GW_test_filtered = Dict([b => GW_test[b] for b in keys(GW)])
  return Omega, GW_test_filtered
end

# check with random test sets
function cc_mobius_test_rand(n_trials::Int64)
  N = free_module(ZZ, 2)
  R, (x, y, u) = polynomial_ring(QQ, [:x, :y, :u])
  for _ in 1:n_trials
    d = Dict{CC, Any}()
    for _ in 1:50
      b = N(rand(Int64, 2) .% 10)
      iszero(b) && continue
      if !haskey(d, b)
        r = rand(Int64, 12)
        d[b] = sum(r .* [one(R), x, y, x*y, x^2, y^2, x^3, x^2*y, x*y^2, u, u^2, u^3])
      end
    end
    cc_mobius_test(d)
  end
end