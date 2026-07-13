using Test
using Oscar
using GKMtools
using CY5

println("== gauge_prediction_2 debug ==")
println("This script only tests code in src/gauge_prediction_2.jl.")
println("It does not call any GW localization routines.")

function _debug_has_odd_u_power(expr, u)
  num = numerator(expr)
  den = denominator(expr)
  @assert den == one(parent(num)) || den isa typeof(one(parent(den)))
  odd_terms = Int[]
  for ev in AbstractAlgebra.exponent_vectors(num)
    upow = ev[end]
    isodd(upow) && push!(odd_terms, upow)
  end
  return sort(unique(odd_terms))
end

@testset "Internal helper sanity" begin
  @test CY5._gp2_part_transpose(Int[]) == Int[]
  @test CY5._gp2_part_transpose([3, 1]) == [2, 1, 1]
  @test CY5._gp2_part_size([3, 1]) == 4
  @test CY5._gp2_part_quad([3, 1]) == 10
  @test CY5._gp2_partitions(0) == [Int[]]
  @test sort(CY5._gp2_partitions(3)) == sort([[3], [2, 1], [1, 1, 1]])
  @test CY5._gp2_mobius(1) == 1
  @test CY5._gp2_mobius(2) == -1
  @test CY5._gp2_mobius(4) == 0
  @test CY5._gp2_mobius(6) == 1
end

@testset "Elementary q-side factors" begin
  S, (rho4, rho5) = polynomial_ring(QQ, ["rho4", "rho5"])
  K = fraction_field(S)

  @test CY5._gp2_tildeZ(Int[], (2, 0), (0, -2), rho4, rho5, K) == K(1)
  @test CY5._gp2_tildeZ([1], (2, 0), (0, -2), rho4, rho5, K) == K(1) // (K(1) - K(rho4)^2)

  r0 = CY5._gp2_R_coeffs(Int[], Int[], 0, rho4, rho5, K)
  @test length(r0) == 1
  @test r0[1] == K(1)

  r1 = CY5._gp2_R_coeffs([1], Int[], 1, rho4, rho5, K)
  @test length(r1) == 2
  @test r1[1] == K(1)
end

@testset "Gauge generator and preimage checks" begin
  G = gkm_5d_gauge(3, 2; equiCY=true)
  Qs, Qts = CY5._gp2_gauge_generators(G, 3)

  @test length(Qs) == 3
  @test length(Qts) == 2
  @test Qs[1] == curve_class(G, "v1", "w1")
  @test Qts[1] == curve_class(G, "v1", "v2")

  dmax, emax = CY5._gp2_preimage_bounds(Qs, Qts, Qs[3])
  @test dmax == (0, 0, 1)
  @test emax == (0, 0)

  dmax2, emax2 = CY5._gp2_preimage_bounds(Qs, Qts, Qs[1])
  @test dmax2[1] >= 1
  @test sum(dmax2) + sum(emax2) > 0
end

@testset "Direct omega_beta_gauge_2 smoke tests" begin
  G = gkm_5d_gauge(3, 2; equiCY=true)
  Q = [curve_class(G, "v$i", "w$i") for i in 1:3]
  Qt = [curve_class(G, "v$i", "v$(i+1)") for i in 1:2]
  S, (rho4, rho5) = polynomial_ring(QQ, ["rho4", "rho5"])
  K = fraction_field(S)

  z_q3 = CY5._gp2_build_Z_h2(G, Q[3], 3, 2, rho4, rho5, K; use_formula_as_written=true)
  zero_key = ntuple(_ -> 0, rank(parent(Q[3])))
  z0 = get(z_q3, zero_key, K(0))

  println()
  println("-- Z-series sanity check before log --")
  println("constant term of Z for Q[3] target = ", z0)
  @test z0 == K(1)

  s_qt = try
    CY5.omega_beta_gauge_2(G, Qt[1], 2)
  catch err
    println("omega_beta_gauge_2(Qt[1]) threw: ", sprint(showerror, err))
    nothing
  end
  s_q3 = try
    CY5.omega_beta_gauge_2(G, Q[3], 2)
  catch err
    println("omega_beta_gauge_2(Q[3]) threw: ", sprint(showerror, err))
    nothing
  end
  s_q1 = try
    CY5.omega_beta_gauge_2(G, Q[1], 1)
  catch err
    println("omega_beta_gauge_2(Q[1]) threw: ", sprint(showerror, err))
    nothing
  end

  @test s_qt !== nothing
  @test s_q3 !== nothing
  @test s_q1 !== nothing

  println()
  println("-- Laurent-series outputs (raw) --")
  println("omega_beta_gauge_2(N=3,m=2, Qt[1], gmax=2) = ", s_qt)
  println("omega_beta_gauge_2(N=3,m=2, Q[3],  gmax=2) = ", s_q3)
  println("omega_beta_gauge_2(N=3,m=2, Q[1],  gmax=1) = ", s_q1)
end

@testset "Prediction wrapper agreement" begin
  G = gkm_5d_gauge(3, 2; equiCY=true)
  Q = [curve_class(G, "v$i", "w$i") for i in 1:3]
  Qt = [curve_class(G, "v$i", "v$(i+1)") for i in 1:2]
  R, (t1, t2, t3, t4, t5, u) = polynomial_ring(QQ, ["t1", "t2", "t3", "t4", "t5", "u"])

  raw = try
    CY5.omega_beta_gauge_2(G, Q[3], 2)
  catch
    nothing
  end
  wrapped = try
    CY5.gkm_5d_gauge_prediction_2(G, [t1, t2, t3, t4, t5], u, Q[3], 2)
  catch err
    println("gkm_5d_gauge_prediction_2(Q[3]) threw: ", sprint(showerror, err))
    nothing
  end

  if raw !== nothing && wrapped !== nothing
    converted = CY5._gp2_series_to_user_ring(raw, [t1, t2, t3, t4, t5], u)
    @test wrapped == converted

    wrapped_qt1 = CY5.gkm_5d_gauge_prediction_2(G, [t1, t2, t3, t4, t5], u, Qt[1], 2)
    wrapped_q1 = CY5.gkm_5d_gauge_prediction_2(G, [t1, t2, t3, t4, t5], u, Q[1], 1)
    odd_q1 = _debug_has_odd_u_power(wrapped_q1, u)

    expected_q3 = (QQ(1)
      + (-QQ(1,24)*t4^2 - QQ(1,24)*t5^2) * u^2
      + (QQ(7,5760)*t4^4 + QQ(1,576)*t4^2*t5^2 + QQ(7,5760)*t5^4) * u^4) //
      (t4*t5*u^2)
    expected_qt1 = (-QQ(2)
      + (-QQ(1,6)*t4^2 - QQ(1,2)*t4*t5 - QQ(1,6)*t5^2) * u^2
      + (QQ(1,360)*t4^4 - QQ(1,72)*t4^2*t5^2 + QQ(1,360)*t5^4) * u^4) //
      (t4*t5*u^2)
    expected_q1 = (QQ(5)
      + (-12*t4 - 12*t5) * u
      + (QQ(403,24)*t4^2 + 34*t4*t5 + QQ(403,24)*t5^2) * u^2) //
      (t4*t5*u^2)

    @test wrapped == expected_q3
    @test wrapped_qt1 == expected_qt1
    @test wrapped_q1 == expected_q1

    println()
    println("-- User-ring outputs --")
    println("gkm_5d_gauge_prediction_2(N=3,m=2, Q[3], gmax=2) = ", wrapped)
    println("gkm_5d_gauge_prediction_2(N=3,m=2, Qt[1], gmax=2) = ", wrapped_qt1)
    println("gkm_5d_gauge_prediction_2(N=3,m=2, Q[1], gmax=1) = ", wrapped_q1)
    println("Odd u-powers appearing in numerator for Q[1]: ", odd_q1)
  else
    @test raw !== nothing
    @test wrapped !== nothing
  end
end

@testset "Inspect historically suspicious classes" begin
  R, (t1, t2, t3, t4, t5, u) = polynomial_ring(QQ, ["t1", "t2", "t3", "t4", "t5", "u"])

  cases = [
    ("N=3 m=2 Q[1]", gkm_5d_gauge(3, 2; equiCY=true), curve_class(gkm_5d_gauge(3, 2; equiCY=true), "v1", "w1"), 1),
    ("N=3 m=2 2Q[2]", gkm_5d_gauge(3, 2; equiCY=true), 2 * curve_class(gkm_5d_gauge(3, 2; equiCY=true), "v2", "w2"), 1),
    ("N=3 m=4 Q[2]", gkm_5d_gauge(3, 4; equiCY=true), curve_class(gkm_5d_gauge(3, 4; equiCY=true), "v2", "w2"), 2),
  ]

  for (label, G, beta, gmax) in cases
    pred = try
      CY5.gkm_5d_gauge_prediction_2(G, [t1, t2, t3, t4, t5], u, beta, gmax)
    catch err
      println()
      println(label, " threw: ", sprint(showerror, err))
      nothing
    end
    pred === nothing && continue
    println()
    println(label, " -> ", pred)
    println("odd u-powers in numerator: ", _debug_has_odd_u_power(pred, u))
  end
end

println()
println("All debug checks in gauge_prediction_2_debug.jl completed.")
