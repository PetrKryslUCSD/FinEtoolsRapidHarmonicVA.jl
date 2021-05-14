# Activate/instantiate the Plotting environment.
# Then execute this file.

include("./define_sim.jl")

using PGFPlotsX
using Plotting
using Plotting: reduced_basis_time, reduced_basis_technique, reduced_basis_style
using Plotting: plot_frf, plot_timing
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA


cdir = sim_directory()
stage = "reduced_basis"

the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("wyd_ritz", "modal"), ("two_stage_wyd_ritz", "modal")]
the_methods = [("free", "modal"), ("two_stage_free", "modal"), ("two_stage_wyd_ritz", "modal")]
the_methods = [("free", "modal"), ("wyd_ritz", "modal"), ("two_stage_free", "modal"), ("two_stage_wyd_ritz", "modal")]

for_nmodes = [50, 100, 200, 400]
plots = []
legends = []
for (reduction_method, harmonic_method) in the_methods

    timings = []
    for nmodes in for_nmodes

        sim = define_sim(; nmodes = (harmonic_method == "direct" ? 0 : nmodes), reduction_method = reduction_method, harmonic_method = harmonic_method)
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
        push!(timings, (nmodes = nmodes, number_of_nodes = number_of_nodes, t))
    end
    
    c = [(t.nmodes, t.t) for t in timings] 
    
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

@pgf ax = SemiLogYAxis(
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
#pgfsave(filename, ax)