using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

using Metis
using FinEtools
using DelimitedFiles
import CoNCMOR: CoNCData, transfmatrix, LegendreBasis, SineCosineBasis 
using FinEtools.MeshExportModule

include("./define_sim.jl")
    
function runme(sim)
    prop = retrieve_json(joinpath(sim_directory(), sim))
    
    fens, fes = make_mesh(prop)
    model = make_model(prop)

    Nc = 64
    C = connectionmatrix(model["femm"], count(fens))
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, Nc; alg = :KWAY)
    mor = CoNCData(fens, partitioning)

    File = joinpath(sim_directory(), "geometry" * ".vtk")
    vtkexportmesh(File,  fes.conn,  fens.xyz, MeshExportModule.VTK.T10)

    for cluster in 1:Nc
        partitioning == cluster
        pnodes = findall(partitioning .== cluster) 
        pPfes = FESetP1(reshape(pnodes, length(pnodes), 1))

        File = joinpath(sim_directory(), "partitioning" * "-cluster-$(cluster).vtk")
        vtkexportmesh(File,  pPfes.conn,  fens.xyz, MeshExportModule.VTK.P1)

    end
end

    sim = define_sim(; mesh_n = 3, nmodes = 50, reduction_method = "two_stage_free")
  runme(sim)