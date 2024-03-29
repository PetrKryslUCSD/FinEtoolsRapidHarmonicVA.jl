module PostUtilities

using JSON
using DelimitedFiles
using Statistics
# using StatsBase
using LinearAlgebra
using Statistics

function fixupext(n, ext)
    f, e = splitext(n)
    return f * ext
end

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


function reduced_basis_technique(reduction_method)
    if reduction_method == "free"
        return "FV"
    elseif reduction_method == "two_stage_free"
        return "2SV"
    elseif reduction_method == "two_stage_free_resid"
        return "2SVR"
    elseif reduction_method == "wyd_ritz"
        return "W"
    elseif reduction_method == "lanczos_ritz"
        return "LNC"
    elseif reduction_method == "two_stage_wyd_ritz"
        return "2SW"
    elseif reduction_method == "two_stage_free_enh"
        return "2SVE"
    else
        return ""
    end
end


function reduced_basis_style(reduction_method)
    if reduction_method == "free"
        return ("red", "diamond")
    elseif reduction_method == "two_stage_free"
        return ("olive", "x")
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
    elseif (reduction_method == "conc_reduced")
        return ("black!40!white", "triangle")
    elseif (reduction_method == "two_stage_free_resid")
        return ("orange", "square")
    elseif reduction_method == "two_stage_free_enh"
        return ("magenta", "triangle")
    else
        return ("", "")
    end
end


function reduced_basis_time(reduction_method, tims)
    if reduction_method == "free"
        return sum([tims[k] for k in ["EV problem"]]) 
    elseif reduction_method == "two_stage_free"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Eigenvector reconstruction"]])
    elseif reduction_method == "two_stage_free_resid"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Additional vectors"]])
    elseif reduction_method == "two_stage_free_enhanced"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Additional vectors"]])
    elseif reduction_method == "wyd_ritz"
        return sum([tims[k] for k in ["Factorize stiffness", "Ritz-vector matrix"]])
    elseif reduction_method == "lanczos_ritz"
        return sum([tims[k] for k in ["Factorize stiffness", "Ritz-vector matrix"]])
    elseif reduction_method == "two_stage_wyd_ritz"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "Factorize stiffness", "Ritz-vector matrix"]])
    elseif reduction_method == "two_stage_free_enh"
        return sum([tims[k] for k in ["Partitioning", "Transformation matrix", "Reduced matrices", "EV problem", "Additional vectors"]])
    else
        return 0.0
    end
end

end # PlotUtilities