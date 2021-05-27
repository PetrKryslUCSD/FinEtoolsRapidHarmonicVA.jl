# Activate/instantiate the Plotting environment.
# Then execute the following file line by line.
using FinEtoolsRapidHarmonicVA

include("./define_sim.jl")

using Plotting
using Plotting: reduced_basis_time
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enhanced", "modal"), ] #  ("lanczos_ritz", "modal")
the_methods = [("none", "direct"), ("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]

for mesh_n in [4, 6, 8]
    for nmodes in [25, ]
        
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


#(reduction_method, harmonic_method) = ("free", "modal")

#for nmodes in [200, ]

#    sims = []
#    for mesh_n in [4, 6, 8, 10, ]
#        sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
#        push!(sims, sim)
#    end

#        #plot_frf_amplitudes(sim_directory(), sims, "frf-m$(mesh_n)-n$(nmodes).pdf")

#        #plot_frf_errors([sim0, sim1, sim2, sim3,  ], "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
#        #plot_frf_errors([sim1,], "frf-m$(mesh_n)-n$(nmodes).pdf")
#    plot_timing_reduced_basis(sim_directory(), sims, "timing-reduced_basis-r$(reduction_method)-n$(nmodes).pdf")
#        #plot_timing([sim1, sim3, sim4, sim5], "harmonic_vibration", "timing-harmonic_vibration-m$(mesh_n)-n$(nmodes).pdf")
#end


#(reduction_method, harmonic_method) = ("free", "modal")

#cdir = sim_directory()
#stage = "reduced_basis"

#timings = []
#for nmodes in [25, 50, 100, 200, ]
#    for mesh_n in [4, 6, 8, 10, ]
#        sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
#        prop = retrieve_json(joinpath(cdir, sim))
#        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
#        results = retrieve_json(j)
#        number_of_nodes = 0
#        t = 0.0
#        if stage in keys(results)
#            sd = results[stage] 
#            tims = sd["timing"]
#            t = reduced_basis_time(prop["reduction_method"], tims)
#            number_of_nodes = sd["number_of_nodes"]
#        end
#        push!(timings, (nmodes = nmodes, mesh_n = mesh_n, number_of_nodes = number_of_nodes, t))
#    end
#end
#@show timings

#l = filter(t -> t.nmodes == 25, timings)
#@show meshes = [t.mesh_n for t in l] 
#c = [(3*t.number_of_nodes, t.t) for t in l] 

#using PGFPlotsX
#@pgf p = PGFPlotsX.Plot(Coordinates(c));
#@pgf ax = Axis(
#    {
#        height = "8cm",
#        width = "9cm",
#        ymin = 0,
#        enlargelimits = "upper",
#        enlarge_x_limits = "true",
#        legend_style = {
#            at = Coordinate(0.5, 0.5),
#            anchor = "south",
#            legend_columns = -1
#        },
#        xlabel = "Degrees of freedom [ND]",
#        ylabel = "Time [s]",
#        nodes_near_coords,
#        nodes_near_coords_align={vertical}
#    },
#    p
#);
#display(ax)
#pgfsave(filename, ax)