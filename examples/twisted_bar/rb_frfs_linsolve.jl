# Activate/instantiate the Plotting environment.
# Then execute this file.
using FinEtoolsRapidHarmonicVA


include("./define_sim.jl")

using Plotting
using Plotting.GnuplotUtilities: clear_terminal, plot_frf_errors, plot_frf_amplitudes, saveas
using FinEtoolsRapidHarmonicVA


the_methods = [ ("two_stage_free_enh", "modal"),]

clear_terminal()

sims = []
sim = define_sim(; mesh_n = 8, nmodes = 0, reduction_method = "none", harmonic_method = "direct")
push!(sims, sim)

for linsolve_method in ["minres"]
    for itmax in [5, 10, 20, 40]
        for mesh_n in [8, ]
            for nmodes in [50, ]
                namebase = "$(linsolve_method)_$(itmax)"
                for m in the_methods
                    reduction_method, harmonic_method = m
                    sim = define_sim(; namebase = namebase, mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, harmonic_method = harmonic_method, linsolve_method = linsolve_method, itmax = itmax)
                    push!(sims, sim)
                end

                title = "mesh_n=$(mesh_n)-nmodes=$(nmodes)" * "minres-itmax=10-time=3.25" 
                plot_frf_errors(sim_directory(), sims, "frf-errors-m$(mesh_n)-n$(nmodes).pdf"; title = title)
                #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

                saveas(title)
            end
        end

    end
end