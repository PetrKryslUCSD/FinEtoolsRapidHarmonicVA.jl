using JSON
using DelimitedFiles

#using PGFPlotsX
# using Random
# using Distributions
using Statistics
# using StatsBase
using LinearAlgebra
using Statistics
using FinEtools
using FinEtoolsRapidHarmonicVA   

function correlation(sim1 = "sim1", sim2 = "sim2")
    prop = retrieve_json(sim1)

        # Load the data for the graph of the FRF
    j = joinpath(prop["resultsdir"], sim1 * "-results" * ".json")
    results = retrieve_json(j)
    hvd = results["harmonic_vibration"] 
        # Unwrap the data
    frequencies = hvd["sweep_frequencies"]
    frf = hvd["frf"]
    mf = frf["file"]
    m1 = retrieve_matrix(mf)
    freal = real.(m1)
    fimag = imag.(m1)

    prop = retrieve_json(sim2)

        # Load the data for the graph of the FRF
    j = joinpath(prop["resultsdir"], sim2 * "-results" * ".json")
    results = retrieve_json(j)
    hvd = results["harmonic_vibration"] 
        # Unwrap the data
    frequencies = hvd["sweep_frequencies"]
    frf = hvd["frf"]
    mf = frf["file"]
    m2 = retrieve_matrix(mf)
    freal = real.(m2)
    fimag = imag.(m2)

    return dot(m1, m2) / norm(m1) / norm(m2)
end

