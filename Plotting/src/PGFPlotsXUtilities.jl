using JSON
using DelimitedFiles

using PGFPlotsX
# using Random
# using Distributions
using Statistics
# using StatsBase
using LinearAlgebra
using Statistics
using FinEtools
using FinEtoolsRapidHarmonicVA   

function plot_timing_reduced_basis(cdir, sim_list = ["sim1"], filename = "plot.pdf")
    stage = "reduced_basis"

    c = []

    for sim in sim_list
        prop = retrieve_json(joinpath(cdir, sim))
        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
        results = retrieve_json(j)
        if stage in keys(results)
            sd = results[stage] 
            tims = sd["timing"]
            t = reduced_basis_time(prop["reduction_method"], tims)
            push!(c, (reduced_basis_technique(prop["reduction_method"]), t))
        end
    end
    @pgf p = PGFPlotsX.Plot(Coordinates(c))
    
    @pgf ax = Axis(
        {
            height = "8cm",
            width = "9cm",
            ybar,
            ymin = 0,
            enlargelimits = "upper",
            enlarge_x_limits = "true",
            legend_style = {
                at = Coordinate(0.5, 0.5),
                anchor = "south",
                legend_columns = -1
            },
            ylabel = "Time [s]",
            symbolic_x_coords=[t[1] for t in c],
            xtick = "data",
            xticklabel_style={
                rotate=45,
                anchor="east"
            },
            nodes_near_coords,
            nodes_near_coords_align={vertical}
        },
        p
    )
    display(ax)
    pgfsave(filename, ax)
    true
end


function plot_times_reduced_basis(list_of_sim_lists = [["sim1",],], filename = "plot.pdf")
    stage = "reduced_basis"

    plots = []

    for sim_list in list_of_sim_lists
        c = []
        s = nothing
        for sim in sim_list
            prop = retrieve_json(sim)
            s = reduced_basis_style(prop["reduction_method"])
            # Load the data for the graph of the FRF
            j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
            results = retrieve_json(j)
            if stage in keys(results)
                sd = results[stage] 
                tims = sd["timing"]
                t = reduced_basis_time(prop["reduction_method"], tims)
                nm = sd["number_of_modes"]
                push!(c, (nm, t))
            end
        end
        #@pgf p = PGFPlotsX.Plot(Coordinates(c))
        
        @pgf p = PGFPlotsX.Plot(
        {
        color = s[1],
        mark = s[2]
        },
        Coordinates(c))
        push!(plots, p)
    end
    
    @pgf ax = Axis(
        {
            height = "8cm",
            width = "9cm",
            ymin = 0,
            enlargelimits = "upper",
            enlarge_x_limits = "true",
            legend_style = {
                at = Coordinate(0.5, 0.5),
                anchor = "south",
                legend_columns = -1
            },
            ylabel = "Time [s]"
        },
        plots...
    )
    display(ax)
    pgfsave(filename, ax)
    true
end

function plot_frf_errors(cdir, sim_list = ["sim1"], filename = "plot.pdf", range = [-Inf, Inf])
    _plot_frf(cdir, sim_list, filename, :errors, range)
end

function plot_frf_amplitudes(cdir, sim_list = ["sim1"], filename = "plot.pdf", range = [-Inf, Inf])
    _plot_frf(cdir, sim_list, filename, :amplitudes, range)
end

function _plot_frf(cdir, sim_list = ["sim1"], filename = "plot.pdf", what = :errors, range = [-Inf, Inf])
    objects = []
    
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
        @pgf p = PGFPlotsX.Plot(
        {
        color = "black",
        line_width  = 0.7
        },
        Coordinates([v for v in  zip(direct_frequencies[range_indexes], direct_ampls[range_indexes])])
        )
        push!(objects, p)
        push!(objects, LegendEntry("Direct"))
        
        direct_frequencies, direct_ampls, range_indexes
    end

    siml = sim_list[2:end]
    if what != :errors
        siml = sim_list
        objects = []
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
        @show s, prop["reduction_method"]
        @pgf p = PGFPlotsX.Plot(
        {
        color = s[1],
        mark = s[2],
        mark_repeat = 15,
        line_width = 0.7
        },
        what == :errors ? 
        Coordinates([v for v in  zip(frequencies[range_indexes], abs.(direct_ampls[range_indexes] .- ampls[range_indexes])./abs.(direct_ampls[range_indexes]))]) : 
        Coordinates([v for v in  zip(frequencies[range_indexes], ampls[range_indexes])])

        )
        push!(objects, p)
        push!(objects, LegendEntry(reduced_basis_technique(prop["reduction_method"])))
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


function plot_frf_errors_direct(sim_list = ["sim1"], filename = "plot.pdf")
    what = :errors
    objects = []
    
    direct_frequencies, direct_ampls, top = let sim = sim_list[1]
        prop = retrieve_json(sim)

        # Load the data for the graph of the FRF
        j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
        mesh_n = prop["mesh_n"]
        results = retrieve_json(j)
        hvd = results["harmonic_vibration"] 
        # Unwrap the data
        direct_frequencies = hvd["sweep_frequencies"]
        frf = hvd["frf"]
        mf = frf["file"]
        m = retrieve_matrix(mf)
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
        prop = retrieve_json(sim)
        mesh_n = prop["mesh_n"]
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

#using PGFPlotsX

## using DataFrames
#import LinearAlgebra: norm
## using Random
## using Distributions
#using Statistics
## using StatsBase

#function makefilename(filename)
#	s = replace(replace(filename, ":"=>"_"), " "=>"_") 
#	if (match(r".*\.pdf$", s) == nothing)
#		s = s * ".pdf"
#	end
#	return s
#end

#function plotconvergence(filename, fst, fs1, fs2; toignore= 0)
#	frange = toignore+1:length(fs1)
#	d1 = fs1[frange]
#	d2 = fs2[frange]
#	dt = fst[frange]
	
#	File = makefilename(filename)
#	df = abs.(d2 - dt) ./ abs.(dt)
#	@pgf pt = Plot({very_thick, color="black", enlargelimits = false}, Table([:x => frange, :y => df]))
#	df = abs.(d2 - d1) ./ abs.(d2)
#	@pgf p1 = Plot({very_thick, color="red", enlargelimits = false}, Table([:x => frange, :y => df]))
#	@pgf ax = Axis({
#	    height = "11cm",
#	    width = "14cm",
#	    xlabel = "Frequency",
#	    ylabel = "Relative difference of the frequencies",
#	    grid="minor",
#	    ymode = "log",
#	    enlargelimits = false
#	}, 
#	pt, p1
#	)
#	display(ax)
#	pgfsave(filename, ax)
#	return ax
#end

#function plotimbalance(filename, imbalanceerror; toignore= 0)
#	frange = toignore+1:length(imbalanceerror)
#	d1 = imbalanceerror[frange]
#	File = makefilename(filename)
#	@pgf pt = Plot({very_thick, color="black", enlargelimits = false}, Table([:x => frange, :y => d1]))
#	@pgf ax = Axis({
#	    height = "11cm",
#	    width = "14cm",
#	    xlabel = "Frequency",
#	    ylabel = "Relative imbalance error",
#	    grid="minor",
#	    enlargelimits = false
#	}, 
#	pt
#	)
#	display(ax)
#	pgfsave(filename, ax)
#	return ax
#end

#function ploteigenvalerror(filename, evt, evs, everr; toignore= 0)
#	frange = toignore+1:length(evt)
#	dt = evt[frange]
#	e1 = everr[frange]
#	d1 = evs[frange]
	
#	File = makefilename(filename)
#	@pgf pp = Plot({very_thick, color="black"}, Table([:x => frange, :y => (d1 + e1) ./ dt]))
#	@pgf pm = Plot({very_thick, color="black"}, Table([:x => frange, :y => (d1 - e1) ./ dt]))
#	@pgf pa = Plot({very_thick, color="blue"}, Table([:x => frange, :y => d1 ./ dt]))
#	@pgf pt = Plot({very_thick, color="red"}, Table([:x => frange, :y => dt ./ dt]))
#	@pgf ax = Axis({
#	    height = "11cm",
#	    width = "14cm",
#	    xlabel = "Frequency",
#	    ylabel = "Normalized eigenvalue",
#	    grid="minor"
#	}, 
#	pp, pm, pt, pa
#	)
#	display(ax)
#	pgfsave(filename, ax)
#	return ax
#end


    
#scoord(key) = begin
#    if key == "Problem setup"
#        return "PS"
#    elseif key == "EV problem"
#        return "EVP"
#    elseif key == "Total"
#        return "Total"       
#    elseif key == "Transformation matrix"
#        return "TM"       
#    elseif key == "Partitioning"
#        return "PART"       
#    elseif key == "Reduced matrices"
#        return "RMK"       
#    elseif key == "Factorize stiffness"
#        return "FSM"       
#    elseif key == "Ritz-vector matrix"
#        return "RVM"       
#    end
#end