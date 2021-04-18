# Activate/instantiate the Plotting environment.
# Then execute the following file line by line.
using FinEtoolsRapidHarmonicVA

include("./define_sim.jl")

using Plotting
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal")]
#methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal")]

for mesh_n in [8, ]
    for nmodes in [25, ]
        
        sims = []
        for (reduction_method, harmonic_method) in methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
            push!(sims, sim)
        end
        
        plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")
        plot_frf_errors(sim_directory(), sims, "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        
        #plot_frf_errors([sim0, sim1, sim2, sim3,  ], "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        #plot_frf_errors([sim1,], "frf-m$(mesh_n)-n$(nmodes).pdf")
        plot_timing_reduced_basis(sim_directory(), sims, "timing-reduced_basis-m$(mesh_n)-n$(nmodes).pdf")
        #plot_timing([sim1, sim3, sim4, sim5], "harmonic_vibration", "timing-harmonic_vibration-m$(mesh_n)-n$(nmodes).pdf")
    end
end