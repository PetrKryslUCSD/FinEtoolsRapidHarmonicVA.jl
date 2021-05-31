module PlotlyJSUtilities

using JSON
using DelimitedFiles
using PlotlyJS
using ..PostUtilities
using ..PostUtilities: reduced_basis_style, reduced_basis_technique
# using Random
# using Distributions
using Statistics
# using StatsBase
using LinearAlgebra
using Statistics
using FinEtools
using FinEtoolsRapidHarmonicVA   

function plot_frf_errors(cdir, sim_list = ["sim1"], filename = "plot.pdf", range = [-Inf, Inf])
    _plot_frf(cdir, sim_list, filename, :errors, range)
end

function plot_frf_amplitudes(cdir, sim_list = ["sim1"], filename = "plot.pdf", range = [-Inf, Inf])
    _plot_frf(cdir, sim_list, filename, :amplitudes, range)
end

colors = Dict("black" => "rgb(15, 15, 15)", "red" => "rgb(160, 15, 15)", "green" => "rgb(15, 160, 15)", "blue" => "rgb(15, 15, 160)", "cyan" => "rgb(15, 160, 160)", "magenta" => "rgb(160, 15, 160)", "orange" => "rgb(160, 130, 60)")

function _plot_frf(cdir, sim_list = ["sim1"], filename = "plot.pdf", what = :errors, range = [-Inf, Inf])
    amplcurves = AbstractTrace[]
    
    direct_frequencies, direct_ampls, range_indexes = let sim = sim_list[1]
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
        range_indexes = [i for i in 1:length(direct_frequencies) 
            if range[1] <= direct_frequencies[i] <= range[2]] 
        c = scatter(;x=direct_frequencies[range_indexes], y=direct_ampls[range_indexes], mode="lines", name = "Reference", line_color = colors["black"], line_width = 2)
        c = scatter(;x=rand(3), y=rand(3), mode="lines", name = "Reference", line_width = 2)
        amplcurves = vcat(amplcurves, c)
        direct_frequencies, direct_ampls, range_indexes
    end

    siml = sim_list[2:end]
    if what != :errors
        siml = sim_list
        amplcurves = AbstractTrace[]
    end
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
        a = reduced_basis_technique(prop["reduction_method"])
        if what == :errors 
            c = scatter(;x=frequencies[range_indexes], y=abs.(direct_ampls[range_indexes] .- ampls[range_indexes]), mode="lines", name = a, line_color = colors[s[1]], line_width = 2)
        else
            c = scatter(;x=frequencies[range_indexes], y=ampls[range_indexes], mode="lines", name = a, line_color = colors[s[1]], line_width = 2)
        end
        amplcurves = vcat(amplcurves, c)
    end

    #layout = Layout(;autosize = true, xaxis=attr(title="Frequency [hertz]", type = "log"), yaxis=attr(title="FRF [mm]", type = "log"))
    layout = Layout(scene = attr(xaxis=attr(title="Frequency [hertz]"), yaxis=attr(title="FRF [mm]")))
    pl = plot(amplcurves, layout; options = Dict(
        :showSendToCloud=>true, 
        :plotlyServerURL=>"https://chart-studio.plotly.com"
        ))

    display(pl)

    #j = joinpath(".", filename)
    #savefig(pl, j)

    true
end


#function plot_frf_errors_direct(sim_list = ["sim1"], filename = "plot.pdf")
#    what = :errors
#    amplcurves = AbstractTrace[]
    
#    direct_frequencies, direct_ampls, top = let sim = sim_list[1]
#        prop = retrieve_json(sim)

#        # Load the data for the graph of the FRF
#        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
#        mesh_n = prop["mesh_n"]
#        results = retrieve_json(j)it
#        hvd = results["harmonic_vibration"] 
#        # Unwrap the data
#        direct_frequencies = hvd["sweep_frequencies"]
#        frf = hvd["frf"]
#        mf = frf["file"]
#        m = retrieve_matrix(mf)
#        freal = real.(m)
#        fimag = imag.(m)
#        # Amplitude graph
#        direct_ampls = abs.(freal + 1im*fimag)./phun("mm")
#        @assert length(direct_ampls) == length(direct_frequencies)
#        @pgf p = PGFPlotsX.Plot(
#        {
#        color = "black"
#        },
#        Coordinates([v for v in  zip(direct_frequencies, direct_ampls)])
#        )
#        push!(objects, p)
#        top = "$(mesh_n)"
#        push!(objects, LegendEntry(top))
        
#        direct_frequencies, direct_ampls, top
#    end

#    for sim in sim_list[2:end]
#        prop = retrieve_json(sim)
#        mesh_n = prop["mesh_n"]
#        # Load the data for the graph of the FRF
#        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
#        results = retrieve_json(j)
#        hvd = results["harmonic_vibration"] 
#        # Unwrap the data
#        frequencies = hvd["sweep_frequencies"]
#        frf = hvd["frf"]
#        mf = frf["file"]
#        m = retrieve_matrix(mf)
#        freal = real.(m)
#        fimag = imag.(m)
#        # Amplitude graph
#        ampls = abs.(freal + 1im*fimag)./phun("mm")
#        @assert length(ampls) == length(frequencies)
#        s = reduced_basis_style(prop["reduction_method"])

#        @pgf p = PGFPlotsX.Plot(
#        {
#        color = "black",
#        text_mark = "$(mesh_n)",
#        mark = "text",
#        mark_repeat = 10+mesh_n
#        },
#        what == :errors ? 
#        Coordinates([v for v in  zip(frequencies, abs.(direct_ampls.-ampls))]) : 
#        Coordinates([v for v in  zip(frequencies, ampls)])
#        )
#        push!(objects, p)
#        leg = top * "-$(mesh_n)"
#        push!(objects, LegendEntry(leg))
#    end

#    @pgf ax = Axis(
#        {
#            xlabel = "Frequency [hertz]",
#            ylabel = "Displacement FRF [mm]",
#            xmode = "log",
#            ymode = "log",
#            yminorgrids = "true",
#            grid = "both",
#            legend_style = {
#                at = Coordinate(0.5, 1.05),
#                anchor = "south",
#                legend_columns = -1
#            },
#        },
#        objects...
#    )

#    display(ax)
#    pgfsave(filename, ax)
#    true
#end



##99999999999999999999
##99999999999999999999
##99999999999999999999
##99999999999999999999
##99999999999999999999

##function plot_frf(sim_list = ["sim1"], filename = "plot.pdf")
##    colors = ["rgb(15, 15, 15)", "rgb(160, 15, 15)", "rgb(15, 160, 15)", "rgb(15, 15, 160)", "rgb(15, 160, 160)"]
##    amplcurves = AbstractTrace[]

##    for sim in sim_list
##        prop = retrieve_json(sim)
##        # Load the data for the graph of the FRF
##        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
##        results = retrieve_json(j)
##        hvd = results["harmonic_vibration"] 
##        # Unwrap the data
##        frequencies = hvd["sweep_frequencies"]
##        frf = hvd["frf"]
##        mf = frf["file"]
##        m = retrieve_matrix(mf)
##        freal = real.(m)
##        fimag = imag.(m)
##        # Amplitude graph
##        ampls = abs.(freal + 1im*fimag)./phun("mm")
##        @assert length(ampls) == length(frequencies)
##        c = scatter(;x=frequencies, y=ampls, mode="lines", name = "$(sim)", line_color = pop!(colors), line_width = 2)
##        push!(amplcurves,  c)
##    end
    
##    # Plot the amplitude of the FRF.

##    # 
##    layout = Layout(;autosize = true, xaxis=attr(title="Frequency [hertz]", type = "log"), yaxis=attr(title="Displacement Amplitude FRF [mm]", type = "log"))
##    # Plot the graphs:
##    #c1 = scatter(;x=frequencies, y=umidAmpl, mode="lines", name = "$(sensortag)", line_color = "rgb(15, 15, 15)")
##    pl = plot(amplcurves, layout; options = Dict(
##            :showSendToCloud=>true, 
##            :plotlyServerURL=>"https://chart-studio.plotly.com"
##            ))
##    display(pl)

##    j = joinpath(".", filename)
##    savefig(pl, j)

    
##    true
##end

##function plot_timing(sim_list = ["sim1"], stage = "harmonic_vibration", filename = "plot.pdf")
##    colors = ["rgb(15, 15, 15)", "rgb(160, 15, 15)", "rgb(15, 160, 15)", "rgb(15, 15, 160)", "rgb(15, 160, 160)"]
##    traces = AbstractTrace[]

##    for sim in sim_list
##        prop = retrieve_json(sim)
##        # Load the data for the graph of the FRF
##        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
##        results = retrieve_json(j)
##        if stage in keys(results)
##            hvd = results[stage] 
##                    # Unwrap the data
##            tims = hvd["timing"]
##            kys = keys(tims)
##            ts = [tims[k] for k in keys(tims)] 
##            c = bar(;x=kys, y=ts, line_color = pop!(colors), name=sim)
##            push!(traces,  c)
##        end
##    end
    
##    # Plot the amplitude of the FRF.

##    # 
##    layout = Layout(;autosize = true, xaxis=attr(title="Operation [ND]"), yaxis=attr(title="Time [s]"))
##    # Plot the graphs:
##    #c1 = scatter(;x=frequencies, y=umidAmpl, mode="lines", name = "$(sensortag)", line_color = "rgb(15, 15, 15)")
##    pl = plot(traces, layout; options = Dict(
##            :showSendToCloud=>true, 
##            :plotlyServerURL=>"https://chart-studio.plotly.com"
##            ))
##    display(pl)
##    j = joinpath(".", filename)
##    savefig(pl, j)
##    true
##end

end # module PlotlyJSUtilities