# Activate/instantiate the Plotting environment.
# Then execute this file.
using FinEtoolsRapidHarmonicVA


include("./define_sim.jl")

using Plotting
using Plotting.PGFPlotsXUtilities: clear_terminal, plot_frf_errors, plot_frf_amplitudes, saveas, plot_frf_errors_direct
# using Plotting.GnuplotUtilities: clear_terminal, plot_frf_errors, plot_frf_amplitudes, saveas# 
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"),]
#the_methods = [("wyd_ritz", "modal"), ("free", "modal"), ("two_stage_free", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal"),]

clear_terminal()

sims = []
for (reduction_method, harmonic_method) in the_methods
    for mesh_n in [8, 6, 4]
        sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
        push!(sims, sim)
    end

    title = "twisted_bar-frf-direct-8-6-4"
    file = title * "-errors" * ".pdf"
    plot_frf_errors_direct(sim_directory(), sims, file)
    file = title * "-amplitudes" * ".pdf"
    plot_frf_errors_direct(sim_directory(), sims, file, what = :amplitudes)
            # plot_frf_amplitudes(sim_directory(), sims, file; title = title)
        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

    saveas(file)
end
