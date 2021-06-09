# Activate/instantiate the Plotting environment.
# Then execute this file.1
using FinEtoolsRapidHarmonicVA


include("./define_sim.jl")

using Plotting
using Plotting.PGFPlotsXUtilities: clear_terminal, plot_frf_errors, plot_frf_amplitudes, saveas
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enhanced", "modal"), ] #  ("lanczos_ritz", "modal")
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]

the_methods = [("free", "modal"), ("two_stage_free", "modal"),]
#the_methods = [("wyd_ritz", "modal"), ("free", "modal"), ("two_stage_free", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal"),]

clear_terminal()



for mesh_n in [3]
    for nmodes in [100, ]
        
        sims = []
        for (reduction_method, harmonic_method) in the_methods
            sim = define_sim(; namebase = "long", mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method,)
            push!(sims, sim)
        end
        
        title = "brake_disc-frf-mesh_n-$(mesh_n)-$(nmodes)-long"
        file = title
        # plot_frf_errors(sim_directory(), sims, file; title = title)
        plot_frf_amplitudes(sim_directory(), sims, file * ".pdf"; range = [1200, 3200], title = title, logx = false)
        # plot_frf_amplitudes(sim_directory(), sims, file * ".pdf"; title = title, logx = false)
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

        # saveas(file  * ".gp")#
        saveas(file * ".pdf")
    end
end


linsolve_method = "minres"
itmax = 20
mesh_n = 3
nmodes = 100
# 
sims = []
(reduction_method, harmonic_method)  = ("free", "modal")
sim = define_sim(; namebase = "long", mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method,)
push!(sims, sim)
(reduction_method, harmonic_method)  = ("two_stage_free_enh", "modal")
sim = define_sim(; namebase = "short", mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
push!(sims, sim)

title = "brake_disc-frf-mesh_n-$(mesh_n)-$(nmodes)-w-enh"
file = title
        # plot_frf_errors(sim_directory(), sims, file; title = title)
plot_frf_amplitudes(sim_directory(), sims, file * ".pdf"; range = [1200, 3200], title = title, logx = false)
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

        # saveas(file  * ".gp")#
saveas(file * ".pdf")

