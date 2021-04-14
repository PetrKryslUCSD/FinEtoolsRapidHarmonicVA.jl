using JSON
using DelimitedFiles
using PlotlyJS
using LinearAlgebra
using Statistics
using FinEtools
using FinEtoolsRapidHarmonicVA   

function plot_frf(sim_list = ["sim1"], filename = "plot.pdf")
    colors = ["rgb(15, 15, 15)", "rgb(160, 15, 15)", "rgb(15, 160, 15)", "rgb(15, 15, 160)", "rgb(15, 160, 160)"]
    amplcurves = AbstractTrace[]

    for sim in sim_list
        prop = retrieve_json(sim)
        # Load the data for the graph of the FRF
        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        m = retrieve_matrix(mf)
        freal = real.(m)
        fimag = imag.(m)
        # Amplitude graph
        ampls = abs.(freal + 1im*fimag)./phun("mm")
        @assert length(ampls) == length(frequencies)
        c = scatter(;x=frequencies, y=ampls, mode="lines", name = "$(sim)", line_color = pop!(colors), line_width = 2)
        push!(amplcurves,  c)
    end
    
    # Plot the amplitude of the FRF.

    # 
    layout = Layout(;autosize = true, xaxis=attr(title="Frequency [hertz]", type = "log"), yaxis=attr(title="Displacement Amplitude FRF [mm]", type = "log"))
    # Plot the graphs:
    #c1 = scatter(;x=frequencies, y=umidAmpl, mode="lines", name = "$(sensortag)", line_color = "rgb(15, 15, 15)")
    pl = plot(amplcurves, layout; options = Dict(
            :showSendToCloud=>true, 
            :plotlyServerURL=>"https://chart-studio.plotly.com"
            ))
    display(pl)

    j = joinpath(".", filename)
    savefig(pl, j)

    
    true
end

function plot_timing(sim_list = ["sim1"], stage = "harmonic_vibration", filename = "plot.pdf")
    colors = ["rgb(15, 15, 15)", "rgb(160, 15, 15)", "rgb(15, 160, 15)", "rgb(15, 15, 160)", "rgb(15, 160, 160)"]
    traces = AbstractTrace[]

    for sim in sim_list
        prop = retrieve_json(sim)
        # Load the data for the graph of the FRF
        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        if stage in keys(results)
            hvd = results[stage] 
                    # Unwrap the data
            tims = hvd["timing"]
            kys = keys(tims)
            ts = [tims[k] for k in keys(tims)] 
            c = bar(;x=kys, y=ts, line_color = pop!(colors), name=sim)
            push!(traces,  c)
        end
    end
    
    # Plot the amplitude of the FRF.

    # 
    layout = Layout(;autosize = true, xaxis=attr(title="Operation [ND]"), yaxis=attr(title="Time [s]"))
    # Plot the graphs:
    #c1 = scatter(;x=frequencies, y=umidAmpl, mode="lines", name = "$(sensortag)", line_color = "rgb(15, 15, 15)")
    pl = plot(traces, layout; options = Dict(
            :showSendToCloud=>true, 
            :plotlyServerURL=>"https://chart-studio.plotly.com"
            ))
    display(pl)
    j = joinpath(".", filename)
    savefig(pl, j)
    true
end