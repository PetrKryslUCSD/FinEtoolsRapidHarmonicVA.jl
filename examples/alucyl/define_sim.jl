using FinEtools
using FinEtoolsRapidHarmonicVA

function alu_cyl_parameters()
    # Material parameters of the solid cylinder
    
    E = 70000*phun("MPa")::FFlt;
    nu = 0.33::FFlt;
    rho = 2700*phun("KG/M^3")::FFlt;
    radius = 0.5*phun("ft"); 
    leng = 2*phun("ft"); 
    mass_shift = (2*pi*100) ^ 2; # to resolve rigid body modes
    
    return E, nu, rho, radius, leng, mass_shift
end

function alu_cyl_mesh(meshmm = 30)
    if meshmm == 30
        return "cylinder-30mm-2116el.mesh"
    elseif meshmm == 20
        return "cylinder-20mm-4436el.mesh"
    elseif meshmm == 15
        return "cylinder-15mm-8204el.mesh"
    elseif meshmm == 10
        return "cylinder-10mm-17980el.mesh"
    else
        return "no such mesh"
    end
end

function alu_cyl_common_parameters()
    prop = Dict() 

    # Folders
    prop["meshesdir"] = "./meshes"
    prop["resultsdir"] = "./results"
    prop["graphicsdir"] = "./graphics"
    prop["matricesdir"] = "./matrices"

    # Material Parameters, geometry
    E, nu, rho, radius, leng, mass_shift = alu_cyl_parameters() 
    
    prop["E"] = E
    prop["nu"] = nu
    prop["rho"] = rho
    prop["mass_shift"] = mass_shift

    # Force and Sensor locations and directions
    prop["force_location"] = [0.0, 0.15227653576, 0.30696353672783677]
    prop["force_direction"] = 2
    
    prop["sensor_location"] = [0.15227653576, 0.0, -0.30696353672783677]
    prop["sensor_direction"] = 2

    # Frequency Sweep
    prop["frequency_sweep"] = (1000, 10000, 1000) # from, to, how many

    return prop
end


function define_sim(; kws...)
    prop = alu_cyl_common_parameters() 
    
    # Defaults
    prop["namebase"] = "sim"
    prop["mesh_mm"] = 15
    prop["neigvs"] = 250
    prop["alpha"] = 1.5
    prop["reduction_method"] = "free_reduced"
    prop["harmonic_method"] = "modal"

    # Overrides
    for k in keys(kws)
        prop[String(k)] = kws[k]
    end

    prop["mesh"] = alu_cyl_mesh(prop["mesh_mm"])
    prop["fmax"] = 2*prop["frequency_sweep"][2]
        
    sim = "$(prop["namebase"])_mesh_$(prop["mesh_mm"])_nev_$(prop["neigvs"])_$(prop["reduction_method"])_$(prop["harmonic_method"])"

    store_json(sim * ".json", prop)

    @info  "Wrote $(sim) json"

    sim
end
