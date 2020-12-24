using JSON
using DelimitedFiles
using PlotlyJS
using LinearAlgebra
using Statistics
using FinEtools
using FinEtools.MeshExportModule


function plot_mesh(sim)
    prop = retrieve_json(sim)
    meshfile = prop["mesh"]
    fens, fes = load_mesh(joinpath(prop["meshesdir"], meshfile))
    vtkfile = with_extension(sim * "-mesh", "vtk")
    mkpath(prop["graphicsdir"])
    vtkexportmesh(joinpath(prop["graphicsdir"], vtkfile), fens, fes)
end

function plot_frf(sim_list = ["sim1"])
    colors = ["rgb(15, 15, 15)", "rgb(125, 15, 15)", "rgb(15, 125, 15)", "rgb(15, 15, 125)"]
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
        c = scatter(;x=frequencies, y=ampls, mode="lines", name = "$(sim)", line_color = pop!(colors))
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
end