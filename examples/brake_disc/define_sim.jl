using LinearAlgebra
using FinEtools
using FinEtools.MeshImportModule
using FinEtoolsDeforLinear
using FinEtoolsRapidHarmonicVA

function sim_directory()
    return dirname(@__FILE__())
end

function make_mesh(prop)
    
    f = joinpath(sim_directory(), "meshes", prop["mesh_name"])
    output = MeshImportModule.import_ABAQUS(f)
    
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm")
    tolerance = 12.5*phun("mm")/1000

    # Clamped face of the shaft hole
    #cpl = output["nsets"]["CIRCULAR-PATCH"]
    clampedl = output["nsets"]["SET-SHAFT"]
    
    # Traction on the circular sub-area of the disk 
    boundaryfes  =   meshboundary(fes);
    loadedl = selectelem(fens, boundaryfes, withnodes = output["nsets"]["SET-CIRCULAR-PATCH"]);
        
    # Sensor Location
    sl = output["nsets"]["SET-CIRCUMFERENCE-EDGE"]
    p = sortperm(fens.xyz[sl, 1])
    sensorn = sl[p][1]

    return fens, fes, clampedl, subset(boundaryfes, loadedl), sensorn
end

function make_model(prop)
    
    fens, fes, clampedl, lbdfes, sensorn = make_mesh(prop)

    # Make fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    setebc!(u, clampedl, true, collect(1:3), 0.0)
    applyebc!(u)
    numberdofs!(u)
    
    # Assuming tetrahedra
    femmmake = FEMMDeforLinear
    integrationrulestiff = TetRule(4)
    integrationrulemass = TetRule(4)

    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    # Compute the stiffness and mass matrices
    femm = femmmake(MR, IntegDomain(fes, integrationrulestiff), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K .= 0.5 * (K .+ transpose(K))

    femm = femmmake(MR, IntegDomain(fes, integrationrulemass), material)
    femm = associategeometry!(femm, geom)
    M = mass(femm, geom, u)
    M .= 0.5 * (M .+ transpose(M))

    # Damping parameters: structural damping assumed
    loss_factor = 0.001
    C = nothing

    fi = ForceIntensity(FFlt[0.0, 0.1*phun("MPa"), 0.1*phun("MPa")]);
    loadedbfemm  = FEMMBase(IntegDomain(lbdfes, TriRule(3)))
    F = distribloads(loadedbfemm, geom, u, fi, 2);
    
    model = Dict()
    model["fens"] = fens
    model["fes"] = fes
    model["femm"] = femm
    model["geom"] = geom
    model["u"] = u
    model["K"] = K
    model["M"] = M
    model["loss_factor"] = loss_factor
    model["C"] = C
    model["F"] = F
    model["sensor_node"] = sensorn
    model["partitioning_method"] = "metis"

    return model
end # 

function exportgraphics(cdir, sim)
    prop = retrieve_json(joinpath(cdir, sim))
    fens, fes, clampedl, lbdfes = make_mesh(prop)
    File =  "twistedbar" * "-nn-$(count(fens))-mesh.vtk"
    vtkexportmesh(joinpath(cdir, File), fens, fes)
    File =  "twistedbar" * "-nn-$(count(fens))-loaded-face-mesh.vtk"
    vtkexportmesh(joinpath(cdir, File), fens, lbdfes)
end

function material_parameters()
    # Material parameters of the structure
    
    E = 210000*phun("MPa")::FFlt;
    nu = 0.3::FFlt;
    rho = 7850*phun("KG/M^3")::FFlt;
    dimensions = [250, 15]*phun("mm"); 
    mass_shift = 0; # to resolve rigid body modes
    
    return E, nu, rho, dimensions, mass_shift
end

function common_parameters()
    prop = Dict() 

    # Folders
    prop["meshesdir"] = "meshes"
    prop["resultsdir"] = "results"
    prop["graphicsdir"] = "graphics"
    prop["matricesdir"] = "matrices"

    # Material Parameters, geometry
    E, nu, rho, dimensions, mass_shift = material_parameters() 
    
    prop["E"] = E
    prop["nu"] = nu
    prop["rho"] = rho
    prop["dimensions"] = dimensions
    prop["mass_shift"] = mass_shift
    prop["smallestdimension"] = minimum(dimensions)
    prop["nbf1maxclamp"] = (3, 6)
    prop["alpha"] = 1.25

    # Sensor location
    #prop["sensor_location"] = [dimensions[1], 0.0, 0.0]
    prop["sensor_direction"] = 3

    # Frequency Sweep
    prop["frequency_sweep"] = (1000, 20000, 1000) # from, to, how many

    return prop
end

function define_sim(; kws...)
    prop = common_parameters() 
    
    # Defaults
    prop["namebase"] = "sim"
    prop["mesh_n"] = 1
    prop["nmodes"] = 60
    prop["reduction_method"] = "free_reduced"
    prop["harmonic_method"] = "modal"
    prop["itmax"] = 0
    prop["resonance_list"] = 1:6
    prop["linsolve_method"] = ""

    # Overrides
    for k in keys(kws)
        prop[String(k)] = kws[k]
    end

    mesh_name = Dict(
        1 => "b250-brake-disk-10096.inp",   
        2 => "b250-brake-disk-21706.inp",
        3 => "b250-brake-disk-52463.inp",
        4 => "b250-brake-disk-130051.inp",  
        5 => "b250-brake-disk-246934.inp",
        )

    prop["mesh_name"]  = mesh_name[prop["mesh_n"]]

    prop["fmax"] = 2*prop["frequency_sweep"][2]

    # Generate the name 
    if prop["reduction_method"] == "two_stage_free_enh"
        @assert prop["linsolve_method"] != "" "Linear system of equation solver must be provided"
        sim = "$(prop["namebase"])_mesh_$(prop["mesh_n"])_nmodes_$(prop["nmodes"])_$(prop["reduction_method"])_$(prop["harmonic_method"])_$(prop["linsolve_method"])_$(prop["itmax"])"
    else # Everyone else
        sim = "$(prop["namebase"])_mesh_$(prop["mesh_n"])_nmodes_$(prop["nmodes"])_$(prop["reduction_method"])_$(prop["harmonic_method"])"
    end

    if (!isfile(sim * ".json")) || 
        ((:force_overwrite in keys(kws)) && kws[:force_overwrite])
        store_json(joinpath(sim_directory(), sim * ".json"), prop)
        @info  "Wrote $(sim) json"
    else
        @warn "Did not overwrite $(sim) json!"
    end
    
    sim
end
