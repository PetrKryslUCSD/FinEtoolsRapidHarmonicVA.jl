module GnuplotUtilities

using JSON
using DelimitedFiles
using Gnuplot
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
    @gp "clear"
    # @gp  :- "set term wxt 0"   :-
    # @gp  :- "set multiplot"  :- 
end

function saveas(a) 
    f, e = splitext(a)
    if e == ".gp"
        Gnuplot.save(a)
    else
        a = fixupext(a, ".pdf")
        @gp  :- "set terminal push" :-
        @gp  :- "set terminal pdfcairo " :-
        @gp  :- "set output '$(a)' "  :-
        @gp  :- "replot"   :-
        @gp  :- "set output" :-
        @gp  :- "set terminal pop" 
    end    
end
    

function plot_frf_errors(cdir, sim_list = ["sim1"], filename = "plot.pdf"; range = [-Inf, Inf], title="")
    _plot_frf(cdir, sim_list, filename, :errors, range, title)
end

function plot_frf_amplitudes(cdir, sim_list = ["sim1"], filename = "plot.pdf"; range = [-Inf, Inf], title="")
    _plot_frf(cdir, sim_list, filename, :amplitudes, range, title)
end

function _plot_frf(cdir, sim_list = ["sim1"], filename = "plot.pdf", what = :errors, range = [-Inf, Inf], title="")
    title = replace(title, '_' => "\\_") 

    @gp  "set terminal windows 0 "  :-
        
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
        range_indexes = [i for i in 1:length(direct_frequencies) if range[1] <= direct_frequencies[i] <= range[2]] 
        plotx = Float64.(direct_frequencies[range_indexes])
        ploty = Float64.(direct_ampls[range_indexes])
        @gp  :- plotx ploty " lw 2 lc rgb 'black' with lines title 'Ref' "  :-
        direct_frequencies, direct_ampls, range_indexes
    end

    siml = sim_list[2:end]
    if what != :errors
        siml = sim_list
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
        plotx = Float64.(direct_frequencies[range_indexes])
        if what == :errors 
            ploty = Float64.(abs.(direct_ampls[range_indexes] .- ampls[range_indexes]))
            @gp  :- plotx ploty " lw 2 lc rgb '$(s[1])' with lines title '$(a)' "  :-
        else
            ploty = Float64.(ampls[range_indexes])
            @gp  :- plotx ploty " lw 2 lc rgb '$(s[1])' with lines title '$(a)' "  :-
        end
    end

    @gp  :- "set logscale x" :-
    @gp  :- "set logscale y" :-
    @gp  :- "set xlabel 'Frequency [Hz]'" :-
    @gp  :- "set ylabel 'FRF [mm]'" :-
    @gp  :- "set title '$(title)'" 


    true
end


end # module GnuplotUtilities