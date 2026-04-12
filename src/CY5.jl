module CY5

using Oscar, GKMtools

import GKMtools.AbstractGKM_graph

include("Types.jl")
include("mobius.jl")
include("strips.jl")
include("closed_vertex.jl")
include("pipeline.jl")
include("brackets.jl")
include("predictions.jl")
include("CY4_vanishing.jl")
include("P1_chain_5d.jl")
include("gauge.jl")
include("gauge_prediction.jl")
include("manual_predictions.jl")
include("manual_examples.jl")

export cc_mobius
export gkm_3d_strip
export gkm_5d_strip
export get_Omega_beta
export gkm_5d_strip_prediction
export Ar_times_C1
export Ar_times_C3
export minus_one_minus_one_chain
export X_times_Ar
export gkm_3d_closed_vertex
export gkm_5d_closed_vertex
export CY5_from_CY4
export gkm_5d_p1_chain
export get_Omega_prediction
export gkm_3d_gauge
export gkm_5d_gauge
export omega_beta_gauge #TODO: remove this export eventually
export gkm_5d_P2_111
export gkm_5d_P3_04
export gkm_5d_P3_13
export gkm_5d_P3_22

end
