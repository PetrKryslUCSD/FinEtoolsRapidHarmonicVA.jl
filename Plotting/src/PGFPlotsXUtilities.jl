module PGFPlotsXUtilities

using JSON
using DelimitedFiles
using PGFPlotsX
using ..PostUtilities
using ..PostUtilities: reduced_basis_style, reduced_basis_technique, fixupext
# using Random
# using Distributions
using Statistics
# using StatsBase
using LinearAlgebra
using Statistics
using FinEtools
using FinEtoolsRapidHarmonicVA   

function clear_terminal()
    # do nothing
end

function saveas(a) 
    f, e = splitext(a)
       
end
    
function plot_frf_errors(cdir, sim_list = ["sim1"], filename = "plot.pdf"; range = [-Inf, Inf], title="")
    _plot_frf(cdir, sim_list, filename, :errors, range, title)
end

function plot_frf_amplitudes(cdir, sim_list = ["sim1"], filename = "plot.pdf"; range = [-Inf, Inf], title="")
    _plot_frf(cdir, sim_list, filename, :amplitudes, range, title)
end

function _plot_frf(cdir, sim_list = ["sim1"], filename = "plot.pdf", what = :errors, range = [-Inf, Inf], title="")
    objects = []
    
    direct_frequencies, direct_ampls = let sim = sim_list[1]
        prop = retrieve_json(joinpath(cdir, sim))

        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        direct_frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        @assert isfile(joinpath(cdir, mf)) "$(joinpath(cdir, mf)) not found"
        m = retrieve_matrix(joinpath(cdir, mf))
        freal = real.(m)
        fimag = imag.(m)
        # Amplitude graph
        direct_ampls = abs.(freal + 1im*fimag)./phun("mm")
        @assert length(direct_ampls) == length(direct_frequencies)
        # range_indexes = [i for i in 1:length(direct_frequencies) 
            # if range[1] <= direct_frequencies[i] <= range[2]] 
        @pgf p = PGFPlotsX.Plot(
        {
        color = "black",
        line_width  = 0.7
        },
        Coordinates([v for v in  zip(direct_frequencies, direct_ampls)])
        )
        push!(objects, p)
        push!(objects, LegendEntry("Ref"))
        
        direct_frequencies, direct_ampls
    end
    if range[1] == -Inf || range[2] == Inf
        range = [minimum(direct_frequencies), maximum(direct_frequencies)]
    end
    siml = sim_list[2:end]
    # if what != :errors
    #     siml = sim_list
    #     # objects = []
    # end
    for sim in siml
        prop = retrieve_json(joinpath(cdir, sim))
        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        m = retrieve_matrix(joinpath(cdir, mf))
        freal = real.(m)
        fimag = imag.(m)
        # Amplitude graph
        ampls = abs.(freal + 1im*fimag)./phun("mm")
        @assert length(ampls) == length(frequencies)
        s = reduced_basis_style(prop["reduction_method"])
        
        @pgf p = PGFPlotsX.Plot(
        {
        color = s[1],
        mark = s[2],
        mark_repeat = 15,
        line_width = 0.7
        },
        what == :errors ? 
        Coordinates([v for v in  zip(frequencies, abs.(direct_ampls .- ampls))]) : 
        Coordinates([v for v in  zip(frequencies, ampls)])

        )
        push!(objects, p)
        push!(objects, LegendEntry(reduced_basis_technique(prop["reduction_method"])))
    end
    
    @pgf ax = Axis(
        {
            xlabel = "Frequency [hertz]",
            ylabel = "Displacement FRF [mm]",
            xmin = range[1],
            xmax = range[2],
            xmode = "log",
            ymode = "log",
            yminorgrids = "true",
            grid = "both",
            legend_style = {
                at = Coordinate(0.5, 1.05),
                anchor = "south",
                legend_columns = -1
            },
        },
        objects...
    )

    display(ax)
    pgfsave(filename, ax)
    true
end


function plot_frf_errors_direct(cdir, sim_list = ["sim1"], filename = "plot.pdf"; what = :errors)
    objects = []
    
    direct_frequencies, direct_ampls, top = let sim = sim_list[1]
        prop = retrieve_json(joinpath(cdir, sim))

        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
        mesh_n = prop["mesh_n"]
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        direct_frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        m = retrieve_matrix(joinpath(cdir, mf))
        freal = real.(m)
        fimag = imag.(m)
        # Amplitude graph
        direct_ampls = abs.(freal + 1im*fimag)./phun("mm")
        @assert length(direct_ampls) == length(direct_frequencies)
        @pgf p = PGFPlotsX.Plot(
        {
        color = "black"
        },
        Coordinates([v for v in  zip(direct_frequencies, direct_ampls)])
        )
        push!(objects, p)
        top = "$(mesh_n)"
        push!(objects, LegendEntry(top))
        
        direct_frequencies, direct_ampls, top
    end

    for sim in sim_list[2:end]
        prop = retrieve_json(joinpath(cdir, sim))
        mesh_n = prop["mesh_n"]
        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        m = retrieve_matrix(joinpath(cdir, mf))
        freal = real.(m)
        fimag = imag.(m)
        # Amplitude graph
        ampls = abs.(freal + 1im*fimag)./phun("mm")
        @assert length(ampls) == length(frequencies)
        s = reduced_basis_style(prop["reduction_method"])

        @pgf p = PGFPlotsX.Plot(
        {
        color = "black",
        text_mark = "$(mesh_n)",
        mark = "text",
        mark_repeat = 10+mesh_n
        },
        what == :errors ? 
        Coordinates([v for v in  zip(frequencies, abs.(direct_ampls.-ampls))]) : 
        Coordinates([v for v in  zip(frequencies, ampls)])
        )
        push!(objects, p)
        leg = top * "-$(mesh_n)"
        push!(objects, LegendEntry(leg))
    end

    @pgf ax = Axis(
        {
            xlabel = "Frequency [hertz]",
            ylabel = "Displacement FRF [mm]",
            xmode = "log",
            ymode = "log",
            yminorgrids = "true",
            grid = "both",
            legend_style = {
                at = Coordinate(0.5, 1.05),
                anchor = "south",
                legend_columns = -1
            },
        },
        objects...
    )

    display(ax)
    pgfsave(filename, ax)
    true
end

end # PGFPlotsXUtilities