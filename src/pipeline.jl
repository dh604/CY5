
@doc raw"""
    get_Omega_beta(G::AbstractGKM_graph, betas::Vector, gMax::Int64; check_predictions::Bool=false)

Return a dictionary assigning to each curve class $\beta$ in `betas` the invariant
$\Omega_b$ as rational function in the equivariant parameters and the genus parameter $u$.

## Arguments
* `G`: the GKM graph of the space $X$ whose $\Omega_\beta$ should be computed.
* `betas`: the vector of curve classes on $X$ for which $\Omega_\beta$ should be computed.
* `gMax`: the maximum genus up to which each $\Omega_\beta$ will be approximated. That is,
          the highest order term will be $u^{2\text{gMax}-2}$.
* `check_predictions`: this optional argument with default value `false` controls whether the computed values of
          $\Omega_\beta$ are compared to the conjectured values. If `check_predictions=true`, the conjectured values
          are retrieved from [`get_Omega_prediction`](@ref).

## Example
TODO
"""
function get_Omega_beta(G::AbstractGKM_graph, betas::Vector{T}, gMax::Int64; check_predictions::Bool=false) where T <: CC
  all_betas = downward_close_ccs(betas)
  res = Dict{CC, Any}()

  # Add u variable for mobius transformation!
  S, t, u = polynomial_ring(QQ, ["t$i" for i in 1:rank_torus(G)], ["u"])
  u = u[1]

  for b in all_betas
    gw0 = gromov_witten(G, b, 0, class_one(); g=0) // 1
    tmp = evaluate(gw0, t)//(u^2)
    for g in 1:gMax
      println("Calculating b=$b, g=$g")
      gw = gromov_witten(G, b, 0, class_one(); g=g) // 1
      tmp += evaluate(gw, t) * u^(2*g - 2)
    end
    res[b] = tmp
  end

  res = cc_mobius(res)

  prediction_tests = Dict{CC, Bool}()

  if check_predictions

    has_CY_substitution = has_attribute(G, :equiCY_substitution)

    if has_CY_substitution
      CY_subst = vcat([evaluate(x, t) for x in get_attribute(G, :equiCY_substitution)], [u])
    end
    
    for b in keys(res)
      p = get_Omega_prediction(G, t, u, b, gMax, res[b])

      if has_CY_substitution
        p = evaluate(p, CY_subst)
      end

      prediction_tests[b] = p == res[b]
      if !prediction_tests[b]
        println("Prediction fails for $b:")
        println("  prediction = $p")
        # println("      actual = $(res[b])")
        get_attribute(G, :example_type) == :gkm_5d_strip_from_3d_CY && println("        type = $(get_attribute(G, :prediction_data)[b])")
      else
        println("Prediction holds for $b")
      end
    end
    if all(b -> prediction_tests[b], keys(prediction_tests))
      println("All predictions hold.")
    else
      println()
    end
  end

  return res
end