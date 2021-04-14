using JSON
using DelimitedFiles
using LinearAlgebra
using Statistics
using FinEtools
using FinEtools.MeshExportModule


function plot_mesh(sim, make_model)
    prop = retrieve_json(sim)
    #meshfile = prop["mesh"]
    #fens, fes = load_mesh(joinpath(prop["meshesdir"], meshfile))
    model = make_model(prop)
    fens, fes = model["fens"], model["fes"]
    vtkfile = with_extension(sim * "-mesh", "vtk")
    mkpath(prop["graphicsdir"])
    vtkexportmesh(joinpath(prop["graphicsdir"], vtkfile), fens, fes)
end
