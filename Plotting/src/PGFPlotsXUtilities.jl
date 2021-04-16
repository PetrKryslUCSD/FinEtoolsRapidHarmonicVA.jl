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

function _reduced_basis_technique(reduction_method)
    if reduction_method == "free"
        return "Full VB"
    elseif reduction_method == "two_stage_free"
        return "2S VB"
    elseif reduction_method == "two_stage_free_enhanced"
        return "Red VB/enh"
    elseif reduction_method == "wyd_ritz"
        return "WYD"
    elseif reduction_method == "lanczos_ritz"
        return "LNC"
    elseif reduction_method == "two_stage_wyd_ritz"
        return "2S WYD"
    else
        return ""
    end
end


function _reduced_basis_style(reduction_method)
    if reduction_method == "free"
        return ("red", "diamond")
    elseif reduction_method == "two_stage_free"
        return ("green", "x")
    elseif reduction_method == "wyd_ritz"
         return ("blue", "o")
    elseif reduction_method == "lanczos_ritz"
         return ("brown", "square")
    elseif reduction_method == "two_stage_free_enhanced"
        return ("magenta", "star")
   elseif reduction_method == "two_stage_wyd_ritz"
        return ("cyan", "star")
    elseif reduction_method == "none"
         return ("black", "+")
    else
        return ("", "")
    end
end

function _reduced_basis_time(reduction_method, tims)
    if reduction_method == "free"
        return sum([tims[k] for k in ["EV problem"]]) 
    elseif reduction_method == "two_stage_free"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Eigenvector reconstruction"]])
    elseif reduction_method == "two_stage_free_enhanced"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Additional vectors"]])
    elseif reduction_method == "wyd_ritz"
        return sum([tims[k] for k in ["Factorize stiffness", "Ritz-vector matrix"]])
    elseif reduction_method == "two_stage_wyd_ritz"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "Factorize stiffness", "Ritz-vector matrix"]])
    else
        return 0.0
    end
end

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
            t = _reduced_basis_time(prop["reduction_method"], tims)
            push!(c, (_reduced_basis_technique(prop["reduction_method"]), t))
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
            @show s = _reduced_basis_style(prop["reduction_method"])
            # Load the data for the graph of the FRF
            j = joinpath(prop["resultsdir"], sim * "-results" * ".json")
            results = retrieve_json(j)
            if stage in keys(results)
                sd = results[stage] 
                tims = sd["timing"]
                t = _reduced_basis_time(prop["reduction_method"], tims)
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

function plot_frf_errors(cdir, sim_list = ["sim1"], filename = "plot.pdf")
    _plot_frf(cdir, sim_list, filename, :errors)
end

function plot_frf_amplitudes(cdir, sim_list = ["sim1"], filename = "plot.pdf")
    _plot_frf(cdir, sim_list, filename, :amplitudes)
end

function _plot_frf(cdir, sim_list = ["sim1"], filename = "plot.pdf", what = :errors)
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
        @show mf
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
        push!(objects, LegendEntry("Direct"))
        
        direct_frequencies, direct_ampls
    end

    for sim in sim_list[2:end]
        prop = retrieve_json(joinpath(cdir, sim))
        # Load the data for the graph of the FRF
        j = joinpath(cdir, prop["resultsdir"], sim * "-results" * ".json")
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
        @show s = _reduced_basis_style(prop["reduction_method"])
        @pgf p = PGFPlotsX.Plot(
        {
        color = s[1],
        mark = s[2],
        mark_repeat = 15,
        },
        what == :errors ? 
        Coordinates([v for v in  zip(frequencies, abs.(direct_ampls.-ampls))]) : 
        Coordinates([v for v in  zip(frequencies, ampls)])

        )
        push!(objects, p)
        push!(objects, LegendEntry(_reduced_basis_technique(prop["reduction_method"])))
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
        s = _reduced_basis_style(prop["reduction_method"])

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