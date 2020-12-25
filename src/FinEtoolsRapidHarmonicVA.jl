module FinEtoolsRapidHarmonicVA 

include("DataUtilities.jl")
export getnamewext, retrieve_json, store_json, retrieve_matrix, store_matrix, with_extension
include("ModelUtilities.jl")
export reduced_basis, show_mesh, harmonic_vibration
include("GraphicsUtilities.jl")
export plot_mesh, plot_frf, plot_timing


end # module
