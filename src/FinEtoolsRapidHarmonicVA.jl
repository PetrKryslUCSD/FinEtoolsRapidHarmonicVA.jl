module FinEtoolsRapidHarmonicVA 

include("DataUtilities.jl")
export getnamewext, retrieve_json, store_json, retrieve_matrix, store_matrix, with_extension
include("ModelUtilities.jl")
export model_setup, full_free_vibration, show_mesh, full_harmonic_vibration, full_harmonic_vibration_direct, redu_free_vibration
include("GraphicsUtilities.jl")
export plot_mesh, plot_frf


end # module
