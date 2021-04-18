# Activate/instantiate the Plotting environment.
# Then execute the following file line by line.

include("./define_sim.jl")

using PGFPlotsX
using Plotting
using Plotting: reduced_basis_time
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


cdir = sim_directory()
stage = "reduced_basis"

the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal")]

plots = []
for (reduction_method, harmonic_method) in the_methods


    timings = []
    for nmodes in [25, 50, 100, 200, ]
        for mesh_n in [4, 6, 8, 10, ]
            sim = define_sim(; mesh_n = mesh_n, nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
            prop = retrieve_json(joinpath(cdir, sim))
            j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
            results = retrieve_json(j)
            number_of_nodes = 0
            t = 0.0
            if stage in keys(results)
                sd = results[stage] 
                tims = sd["timing"]
                t = reduced_basis_time(prop["reduction_method"], tims)
                number_of_nodes = sd["number_of_nodes"]
            end
            push!(timings, (nmodes = nmodes, mesh_n = mesh_n, number_of_nodes = number_of_nodes, t))
        end
    end
    @show timings

    l = filter(t -> t.nmodes == 25, timings)
    meshes = [t.mesh_n for t in l] 
    c = [(3*t.number_of_nodes, t.t) for t in l] 

    @pgf p = PGFPlotsX.Plot(Coordinates(c));
    push!(plots, p)
end

@pgf ax = Axis(
    {
        height = "8cm",
        width = "9cm",
        ymode = "log",
        enlargelimits = "upper",
        enlarge_x_limits = "true",
        legend_style = {
            at = Coordinate(0.5, 0.5),
            anchor = "south",
            legend_columns = -1
        },
        xlabel = "Degrees of freedom [ND]",
        ylabel = "Time [s]",
        nodes_near_coords,
        nodes_near_coords_align={vertical}
    },
    plots...
);
display(ax)
#pgfsave(filename, ax)