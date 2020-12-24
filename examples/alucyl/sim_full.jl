using FinEtools
using FinEtoolsRapidHarmonicVA

include("alu_cyl_parameters.jl")

sim = let
    prop = Dict() 

    # Mesh folder
    prop["meshesdir"] = "./meshes"

    # Output folder
    prop["resultsdir"] = "./results"

    # Graphics folder
    prop["graphicsdir"] = "./graphics"

    # Matrix folder
    prop["matricesdir"] = "./matrices"

    prop["mesh"] = "cylinder-30mm-2116el.mesh"
    #prop["mesh"] = "cylinder-15mm-8204el.mesh"
    #prop["mesh"] = "cylinder-10mm-17980el.mesh"

    E, nu, rho, radius, leng, mass_shift = alu_cyl_parameters() 
    
    prop["E"] = E
    prop["nu"] = nu
    prop["rho"] = rho
    prop["mass_shift"] = mass_shift

    prop["neigvs"] = 150

    prop["force_location"] = [0.0, 0.15227653576, 0.30696353672783677]
    prop["force_direction"] = 2
    
    prop["sensor_location"] = [0.15227653576, 0.0, -0.30696353672783677]
    prop["sensor_direction"] = 2

    prop["frequency_sweep"] = (1000, 10000, 1000) # from, to, how many

    sim = "$(splitext(basename(@__FILE__()))[1])_nev_$(prop["neigvs"])"

    store_json(sim * ".json", prop)

    @info  "Wrote $(sim) json"

    sim
end
