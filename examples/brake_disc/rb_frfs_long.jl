# Activate/instantiate the Plotting environment.
# Then execute this file.1
using FinEtoolsRapidHarmonicVA


include("./define_sim.jl")

using Plotting
using Plotting.PGFPlotsXUtilities: clear_terminal, plot_frf_errors, plot_frf_amplitudes, saveas
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enhanced", "modal"), ] #  ("lanczos_ritz", "modal")
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]
#the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_residual", "modal"),]
the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal"),]
the_methods = [("free", "modal"), ("two_stage_free", "modal"),]
the_methods = [("free", "modal"), ("two_stage_free_enh", "modal"),]

clear_terminal()

linsolve_method = "minres"
itmax = 20
for mesh_n in [3]
    for nmodes in [100, ]
        
        sims = []
        for (reduction_method, harmonic_method) in the_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method, linsolve_method = linsolve_method, itmax = itmax)
            push!(sims, sim)
        end
        
        title = "brake_disc-frf-mesh_n-$(mesh_n)-$(nmodes)-$(linsolve_method)-$(itmax)"
        file = title
        # plot_frf_errors(sim_directory(), sims, file; title = title)
        plot_frf_amplitudes(sim_directory(), sims, file * ".pdf"; title = title, logx = false, range = [1200, 3200], mark_repeat = 159)
        # 
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

        # saveas(file  * ".gp")#
        saveas(file * ".pdf")
    end
end
