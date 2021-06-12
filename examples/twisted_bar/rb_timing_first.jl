# Activate/instantiate the Plotting environment.
# Then run this file.

include("./define_sim.jl")

using PGFPlotsX
using Plotting
using Plotting.PostUtilities: reduced_basis_style, reduced_basis_technique, reduced_basis_time, fixupext
using FinEtoolsRapidHarmonicVA

cdir = sim_directory()
stage = "reduced_basis"

the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ("two_stage_free_enh", "modal")]
the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal"), ]
#the_methods = [("free", "modal"), ("two_stage_free", "modal"), ]

for_nmodes = 25
plots = []
legends = []

for (reduction_method, harmonic_method) in the_methods


    timings = []
    for nmodes in [for_nmodes,  ]
        for mesh_n in [4, 6, 8, 10, 12, 14, 16]
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
    #@show timings

    l = filter(t -> t.nmodes == for_nmodes, timings)
    meshes = [t.mesh_n for t in l] 
    c = [(3*t.number_of_nodes, t.t) for t in l] 

    s = reduced_basis_style(reduction_method)
    @pgf p = PGFPlotsX.Plot({
        color = s[1],
        mark = s[2], 
        mark_size = 4,
        line_width = 2
        }, Coordinates(c));
    push!(plots, p)
    @pgf leg = LegendEntry(reduced_basis_technique(reduction_method))
    push!(legends, leg)
end

@pgf ax = LogLogAxis(
    {
        height = "8cm",
        width = "9cm",
        log_basis_y = 10,
        enlargelimits = "upper",
        enlarge_x_limits = "true",
        legend_style = {
            at = Coordinate(0.5, 1.01),
            anchor = "south",
            legend_columns = -1
        },
        xlabel = "Degrees of freedom [ND]",
        ylabel = "Time [s]",
        nodes_near_coords,
        nodes_near_coords_align={horizontal}
    },
    plots..., legends...
);
display(ax)
pgfsave("twisted_bar-timing-$(for_nmodes)-$(linsolve_method)-$(itmax).pdf", ax)