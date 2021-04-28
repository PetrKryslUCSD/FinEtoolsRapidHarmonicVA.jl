using Metis
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsRapidEig
using Arpack
using Statistics: mean
using FinEtoolsRapidEig.DataUtilities: getjsonname, savejson, getmdatname, store
using LinearAlgebra: issymmetric, norm, ldlt
using SparseArrays
using DelimitedFiles

include("lug_utilities.jl")


function runme(meshfile, neigvs)

	@info "FULL model, mesh $(meshfile)"

	problem, E, nu, rho = lug_parameters()

	properties = Dict("problem"=>problem, "E"=>E, "nu"=>nu, "rho"=>rho, "meshfile"=>meshfile, "neigvs"=>neigvs)
	features = "$(meshfile)-neigvs=$(neigvs)"

	timing = Dict{String, FFlt}()

	timing["Problem setup"] = @elapsed begin
		fens, fesets, integrationrulestiff, integrationrulemass = lug_load_mesh(meshfile)
        u, K, M = lug_setup(E, nu, rho, fens, fesets, integrationrulestiff, integrationrulemass, :consistent)
		# @show issymmetric(K) && issymmetric(M)
    end

	#timing["EV problem"] = @elapsed begin
	#	eval, evec = eigs(K, M; nev=neigvs, maxiter = 600, which=:SM)
	#	fs = real(sqrt.(complex(eval)))/(2*pi)
	#end
	#println("Reference Eigenvalues: $fs [Hz]")

	#properties["frequencies"] = fs
	#timing["Total"] = timing["Problem setup"] + timing["EV problem"]
	#properties["timing"] = timing

	#properties["eigenvectors"] = Dict("dims"=>size(evec), "type"=>"$(eltype(evec))", "file"=>getmdatname(problem, features, "evec"))
	#store(evec, properties["eigenvectors"]["file"])

	#savejson(getjsonname(problem, features), properties)

    fes = fesets[1]
    for i in 2:length(fesets)
        fes = cat(fes, fesets[i])
    end

    C = connectionmatrix(FEMMBase(IntegDomain(fes, TetRule(1))), count(fens))
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, 8; alg = :KWAY)

    @show realNc = maximum(unique(partitioning))
    for cluster = 1:realNc
        pnodes = findall(partitioning .== cluster) 
        @show length(pnodes)
        pPfes = FESetP1(reshape(pnodes, length(pnodes), 1))
        
        File =  meshfile * "-cluster-$(cluster).vtk"
        vtkexportmesh(File,  pPfes.conn,  fens.xyz, FinEtools.MeshExportModule.VTK.P1)
    end

	timing
end

# meshes = [   
# "lug-217666.nas",
# "lug-121531.nas",  
# "lug-54776.nas",   
# "lug-18540.nas"
# ]
# for meshfile in meshes
# 	for neigvs in [15, 50, 150, 300, 600]#
# 		@show meshfile, neigvs
# 		@show timing = runme(meshfile, neigvs)
# 	end
# end

meshes = [   
# "lug-217666.nas",
# "lug-121531.nas",  
# "lug-54776.nas",   
"lug-18540.nas"
]
for meshfile in meshes
	for neigvs in [150]#
		@show meshfile, neigvs
		@show timing = runme(meshfile, neigvs)
	end
end
