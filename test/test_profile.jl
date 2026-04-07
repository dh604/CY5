using Oscar
using CY5

R, (T0, T1, T2, T3, T4, u) = polynomial_ring(QQ, ["T0", "T1", "T2", "T3", "T4", "u"])

println("=== Full omega_series for P2_111, d_max=4, deg=6 ===")
@time result = CY5.omega_series("Omega_P2_111", 4, [T0, T1, T2, T3, T4], u, 6)
for (d, r) in enumerate(result)
  s = string(r)
  println("  d=$d: $(length(s)) chars")
end
