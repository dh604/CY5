using Test
using Oscar
using GKMtools
using GW_CY5

# ==============================================================================
# Sanity check for the gauge prediction in src/gauge_prediction_3.jl.
#
# We compute the *predicted* Omega_beta (via gkm_5d_gauge_prediction_3) over the
# conjecturally valid range |m| < N for a representative set of curve classes,
# and verify two structural properties of every prediction:
#
#   * only EVEN powers of u appear, and
#   * no power of u below -2 appears.
#
# (Omega_beta is a genus expansion 1/u^2 + ... whose nonzero terms sit at the
# even powers u^{2g-2}, g >= 0.)  No actual GW invariants are computed here.
# ==============================================================================

println("== gauge_prediction_3 sanity check (no GW localization) ==")

# The true u-exponents of a prediction, returned as a sorted Vector{Int}.
#
# A prediction lives in the fraction field QQ(t1,...,t5,u); concretely it equals
# num(t,u) / den(t,u) with den a pure power of u times a u-free polynomial.  The
# u-exponents of the value are therefore the u-degrees of the numerator's terms,
# each shifted down by the (single) u-degree of the denominator.  The zero
# prediction has no exponents.
function _gp3_u_exponents(f)
  iszero(f) && return Int[]
  num = numerator(f)
  den = denominator(f)
  den_u = unique(ev[end] for ev in AbstractAlgebra.exponent_vectors(den))
  @assert length(den_u) == 1 "denominator is not a pure power of u: $den"
  d = den_u[1]
  return sort(unique(ev[end] - d for ev in AbstractAlgebra.exponent_vectors(num)))
end

@testset "gp3 prediction u-exponents (|m| < N)" begin
  R, (t1, t2, t3, t4, t5, u) = polynomial_ring(QQ, ["t1", "t2", "t3", "t4", "t5", "u"])
  t = [t1, t2, t3, t4, t5]

  n_checked = 0
  violations = Tuple{Int,Int,Any,Int,Vector{Int}}[]   # (N, m, beta, gMax, bad exps)

  for N in 2:4, m in -(N-1):(N-1)           # the whole range |m| < N
    G = gkm_5d_gauge(N, m; equiCY=true)
    Q  = [curve_class(G, "v$i", "w$i") for i in 1:N]
    Qt = [curve_class(G, "v$i", "v$(i+1)") for i in 1:(N-1)]

    # Always test every section / fibre generator.  For the cheaper small ranks
    # we additionally test a fibre multiple and a mixed sum, which exercise the
    # multi-cover / Moebius path and a non-generator class.
    betas = typeof(Q[1])[]
    append!(betas, Q)
    append!(betas, Qt)
    if N <= 3 && !isempty(Qt)
      push!(betas, 2 * Qt[1])
      push!(betas, Q[1] + Qt[1])
    end
    gMaxes = N <= 3 ? (1:2) : (2:2)

    for beta in betas, gMax in gMaxes
      pred = GW_CY5.gkm_5d_gauge_prediction_3(G, t, u, beta, gMax)
      exps = _gp3_u_exponents(pred)
      bad = filter(e -> isodd(e) || e < -2, exps)
      isempty(bad) || push!(violations, (N, m, beta, gMax, bad))
      n_checked += 1
      @test isempty(bad)
    end
    println("  done N=$N, m=$m  ($(length(betas)) classes)"); flush(stdout)
  end

  println("Checked $n_checked predictions across N in 2:4 and all |m| < N.")
  if isempty(violations)
    println("All predictions use only even powers of u and none below u^-2.")
  else
    println("Sanity VIOLATIONS:")
    for (N, m, beta, gMax, bad) in violations
      println("  N=$N m=$m beta=$beta gMax=$gMax -> offending u-exponents $bad")
    end
  end
  @test isempty(violations)
end

println("gauge_prediction_3_sanity.jl completed.")
