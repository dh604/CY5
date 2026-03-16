# Spaces

On this page, we list the implemented spaces and constructor functions that are compatible with the [Pipeline](pipeline.md).
This means that they can be used in comjunction with [`get_Omega_prediction`](@ref) and hence with [`get_Omega_beta`](@ref) and its optional argument `check_predictions=true`.

!!! note
    Only those functions that return a 5-dimensional result are compatible with the [Pipeline](pipeline.md).
    Functions returning a lower-dimensional result are listed because they are ingredients for 5-dimensional constructions.

## Strip geometries
```@docs
Ar_times_C1
Ar_times_C3
gkm_3d_strip
gkm_5d_strip
minus_one_minus_one_chain
```

## Free $\mathbb{P}^1$ chains
```@docs
gkm_5d_p1_chain
```

## Closed vertex

```@docs
gkm_3d_closed_vertex
gkm_5d_closed_vertex
```

## GW-Gauge correspondence
```@docs
gkm_3d_gauge
gkm_5d_gauge
```

## From 3-folds
```@docs
X_times_Ar
```

## From 4-folds
```@docs
CY5_from_CY4
```