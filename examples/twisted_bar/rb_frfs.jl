# Activate/instantiate the Plotting environment.
# Then execute this file.
using FinEtoolsRapidHarmonicVA

include("./define_sim.jl")

using Plotting
using Plotting: plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enhanced", "modal"), ] #  ("lanczos_ritz", "modal")
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]
#the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_residual", "modal"),]
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal"),]
#the_methods = [("wyd_ritz", "smodal"), ("free", "modal"), ("two_stage_free", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal"),]

for mesh_n in [8]
    for nmodes in [25, ]
        
        sims = []
        for (reduction_method, harmonic_method) in the_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
            push!(sims, sim)
        end
        
        plot_frf_errors(sim_directory(), sims, "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")
    end
end
