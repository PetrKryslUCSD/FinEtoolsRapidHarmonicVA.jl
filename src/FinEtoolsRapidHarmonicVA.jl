module FinEtoolsRapidHarmonicVA 

include("DataUtilities.jl")
export getnamewext, retrieve_json, store_json, retrieve_matrix, store_matrix, with_extension

#include("AssembleReduced.jl")
#export SysmatAssemblerReduced, startassembly!, assemble!, makematrix!

include("ModelUtilities.jl")
export reduced_basis, show_mesh, harmonic_vibration, solve

include("GraphicsUtilities.jl")
export plot_mesh

end # module
