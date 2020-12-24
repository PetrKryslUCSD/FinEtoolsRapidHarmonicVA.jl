module PGFPlotsXUtilities

using PGFPlotsX

# using DataFrames
import LinearAlgebra: norm
# using Random
# using Distributions
using Statistics
# using StatsBase

function makefilename(filename)
	s = replace(replace(filename, ":"=>"_"), " "=>"_") 
	if (match(r".*\.pdf$", s) == nothing)
		s = s * ".pdf"
	end
	return s
end

function plotconvergence(filename, fst, fs1, fs2; toignore= 0)
	frange = toignore+1:length(fs1)
	d1 = fs1[frange]
	d2 = fs2[frange]
	dt = fst[frange]
	
	File = makefilename(filename)
	df = abs.(d2 - dt) ./ abs.(dt)
	@pgf pt = Plot({very_thick, color="black", enlargelimits = false}, Table([:x => frange, :y => df]))
	df = abs.(d2 - d1) ./ abs.(d2)
	@pgf p1 = Plot({very_thick, color="red", enlargelimits = false}, Table([:x => frange, :y => df]))
	@pgf ax = Axis({
	    height = "11cm",
	    width = "14cm",
	    xlabel = "Frequency",
	    ylabel = "Relative difference of the frequencies",
	    grid="major",
	    ymode = "log",
	    enlargelimits = false
	}, 
	pt, p1
	)
	display(ax)
	pgfsave(filename, ax)
	return ax
end

function plotimbalance(filename, imbalanceerror; toignore= 0)
	frange = toignore+1:length(imbalanceerror)
	d1 = imbalanceerror[frange]
	File = makefilename(filename)
	@pgf pt = Plot({very_thick, color="black", enlargelimits = false}, Table([:x => frange, :y => d1]))
	@pgf ax = Axis({
	    height = "11cm",
	    width = "14cm",
	    xlabel = "Frequency",
	    ylabel = "Relative imbalance error",
	    grid="major",
	    enlargelimits = false
	}, 
	pt
	)
	display(ax)
	pgfsave(filename, ax)
	return ax
end

function ploteigenvalerror(filename, evt, evs, everr; toignore= 0)
	frange = toignore+1:length(evt)
	dt = evt[frange]
	e1 = everr[frange]
	d1 = evs[frange]
	
	File = makefilename(filename)
	@pgf pp = Plot({very_thick, color="black"}, Table([:x => frange, :y => (d1 + e1) ./ dt]))
	@pgf pm = Plot({very_thick, color="black"}, Table([:x => frange, :y => (d1 - e1) ./ dt]))
	@pgf pa = Plot({very_thick, color="blue"}, Table([:x => frange, :y => d1 ./ dt]))
	@pgf pt = Plot({very_thick, color="red"}, Table([:x => frange, :y => dt ./ dt]))
	@pgf ax = Axis({
	    height = "11cm",
	    width = "14cm",
	    xlabel = "Frequency",
	    ylabel = "Normalized eigenvalue",
	    grid="major"
	}, 
	pp, pm, pt, pa
	)
	display(ax)
	pgfsave(filename, ax)
	return ax
end

end

# 	    ymode = "log",