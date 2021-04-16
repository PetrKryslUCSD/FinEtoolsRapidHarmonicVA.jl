using FinEtoolsRapidHarmonicVA

include("../examples/twisted_bar/define_sim.jl")

using Plotting
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
reduction_methods = ["free", "two_stage_free"]

for mesh_n in [4, ]
    for nmodes in [50, ]
     
     sims = []
     for reduction_method in reduction_methods
                 sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
                 push!(sims, sim)
             end
        
        #sim5 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced_enhanced", harmonic_method = "modal")
        #sim6 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "reduced_wyd_ritz", harmonic_method = "modal")
        plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")
        #plot_frf_errors([sim0, sim1, sim2, sim3,  ], "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        #plot_frf_errors([sim1,], "frf-m$(mesh_n)-n$(nmodes).pdf")
        #plot_timing_reduced_basis([sim1, sim2, sim3, sim4], "timing-reduced_basis-m$(mesh_n)-n$(nmodes).pdf")
        #plot_timing([sim1, sim3, sim4, sim5], "harmonic_vibration", "timing-harmonic_vibration-m$(mesh_n)-n$(nmodes).pdf")
    end
end