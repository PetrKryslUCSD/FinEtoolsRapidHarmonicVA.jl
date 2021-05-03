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
    clampedl = output["nsets"]["SHAFT"]
    
    # Traction on the circular sub-area of the disk 
    boundaryfes  =   meshboundary(fes);
    loadedl = selectelem(fens, boundaryfes, withnodes = output["nsets"]["CIRCULAR-PATCH"]);
        
    # Sensor Location
    sl = output["nsets"]["CIRCUMFERENCE-EDGE"]
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

    # Compute the parameters of Rayleigh damping. For the two selected
    # frequencies we have the relationship between the damping ratio and
    # the Rayleigh parameters
    a0(zeta1, zeta2, o1, o2) = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);# a0
    a1(zeta1, zeta2, o1, o2) = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);# a1

    # Damping parameters
    zeta1, zeta2, omega1, omega2 = 0.01, 0.005, 100, 5000
    a0v, a1v = a0(zeta1, zeta2, omega1, omega2), a1(zeta1, zeta2, omega1, omega2) 
    # The damping is evaluated on the tetrahedral mesh using the full integration rule
    femmc = FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
    C = a0v * mass(femmc, geom, u) + a1v * stiffness(femmc, geom, u)

    fi = ForceIntensity(FFlt[0.0, 0.0, 0.01*phun("MPa")]);
    loadedbfemm  = FEMMBase(IntegDomain(lbdfes, TriRule(1)))
    F = distribloads(loadedbfemm, geom, u, fi, 2);
    
    model = Dict()
    model["fens"] = fens
    model["fes"] = fes
    model["femm"] = femm
    model["geom"] = geom
    model["u"] = u
    model["K"] = K
    model["M"] = M
    model["C"] = C
    model["F"] = F
    model["sensor_node"] = sensorn

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
    
    E = 70000*phun("MPa")::FFlt;
    nu = 0.33::FFlt;
    rho = 2700*phun("KG/M^3")::FFlt;
    dimensions = [250, 15, 25]*phun("mm"); 
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

    prop["partitioning_method"] = "metis"

    # Material Parameters, geometry
    E, nu, rho, dimensions, mass_shift = material_parameters() 
    
    prop["E"] = E
    prop["nu"] = nu
    prop["rho"] = rho
    prop["dimensions"] = dimensions
    prop["mass_shift"] = mass_shift
    prop["smallestdimension"] = minimum(dimensions)
    prop["nbf1maxclamp"] = (3, 3)
    prop["alpha"] = 1.5

    # Sensor location
    #prop["sensor_location"] = [dimensions[1], 0.0, 0.0]
    prop["sensor_direction"] = 2

    # Frequency Sweep
    prop["frequency_sweep"] = (100, 20000, 200) # from, to, how many

    return prop
end

function define_sim(; kws...)
    prop = common_parameters() 
    
    # Defaults
    prop["namebase"] = "sim"
    prop["mesh_name"] = "b250-brake-disc.inp"
    prop["nmodes"] = 60
    prop["reduction_method"] = "free_reduced"
    prop["harmonic_method"] = "modal"

    # Overrides
    for k in keys(kws)
        prop[String(k)] = kws[k]
    end

    prop["fmax"] = 2*prop["frequency_sweep"][2]
        
    sim = "$(prop["namebase"])_mesh_$(prop["mesh_n"])_nmodes_$(prop["nmodes"])_$(prop["reduction_method"])_$(prop["harmonic_method"])"

    if (!isfile(sim * ".json")) || 
        ((:force_overwrite in keys(kws)) && kws[:force_overwrite])
        store_json(joinpath(sim_directory(), sim * ".json"), prop)
        @info  "Wrote $(sim) json"
    else
        @warn "Did not overwrite $(sim) json!"
    end
    
    sim
end
