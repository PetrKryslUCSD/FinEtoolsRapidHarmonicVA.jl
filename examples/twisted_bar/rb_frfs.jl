# Activate/instantiate the Plotting environment.
# Then execute this file.
using FinEtoolsRapidHarmonicVA

include("./define_sim.jl")

using Plotting
using Plotting: reduced_basis_time
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enhanced", "modal"), ] #  ("lanczos_ritz", "modal")
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]
#the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_residual", "modal"),]
for mesh_n in [8]
    for nmodes in [200, ]
        
        sims = []
        for (reduction_method, harmonic_method) in the_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
            push!(sims, sim)
        end
        
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf", [190, 210])
        #
        plot_frf_errors(sim_directory(), sims, "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")
    end
end


#using Plotting: plot_frf_errors_direct

#    sims = []
#for mesh_n in [8, 6, 4, ]
#    sim = define_sim(; mesh_n = mesh_n, nmodes = 0, reduction_method = "none", harmonic_method = "direct")
#    push!(sims, sim)
#end
#plot_frf_errors_direct(sims, "plot.pdf")