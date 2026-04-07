using Oscar
using CY5

println("=== Test 1: load_omega_data with specific degrees ===")
Qs = CY5.load_omega_data("Omega_P2_111", [1, 3])
println("Keys: ", sort(collect(keys(Qs))))
for d in sort(collect(keys(Qs)))
  s = string(Qs[d])
  println("  Q$d: $(length(s)) chars")
end

println("\n=== Test 2: load_omega_data with d_max (convenience) ===")
Qs2 = CY5.load_omega_data("Omega_P2_111", 3)
println("Keys: ", sort(collect(keys(Qs2))))

println("\n=== Test 3: omega_exp_substitution ===")
Qs1 = CY5.load_omega_data("Omega_P2_111", [1])
polys, t = CY5.omega_exp_substitution(Qs1, 6)
println("Q1 deg=6: ", polys[1])

println("\n=== Test 4: omega_u_scaling ===")
R, (T0, T1, T2, T3, T4, u) = polynomial_ring(QQ, ["T0", "T1", "T2", "T3", "T4", "u"])
scaled = CY5.omega_u_scaling(polys, [T0, T1, T2, T3, T4], u)
println("Keys: ", sort(collect(keys(scaled))))

println("\n=== Test 5: omega_series with collection ===")
result = CY5.omega_series("Omega_P2_111", [1, 3], [T0, T1, T2, T3, T4], u, 6)
println("Keys: ", sort(collect(keys(result))))
for d in sort(collect(keys(result)))
  println("  d=$d: $(length(string(result[d]))) chars")
end

println("\n=== Test 6: omega_series with d_max ===")
result2 = CY5.omega_series("Omega_P2_111", 3, [T0, T1, T2, T3, T4], u, 6)
println("Keys: ", sort(collect(keys(result2))))

println("\n=== Test 7: consistency check ===")
Qs_all = CY5.load_omega_data("Omega_P2_111", 3)
polys_all, _ = CY5.omega_exp_substitution(Qs_all, 6)
scaled_all = CY5.omega_u_scaling(polys_all, [T0, T1, T2, T3, T4], u)
println("d=1 matches: ", scaled_all[1] == result[1])
println("d=3 matches: ", scaled_all[3] == result[3])
