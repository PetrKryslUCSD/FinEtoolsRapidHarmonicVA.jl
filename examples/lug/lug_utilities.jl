function lug_parameters()
	# Material parameters of the load-bearing lug
	problem = "lug"
	E = 200000*phun("MPa")::FFlt;
	nu = 0.33::FFlt;
	rho = 7850*phun("KG/M^3")::FFlt;
	return problem, E, nu, rho
end

function lug_load_mesh(meshfile)
	meshfilebase, ext = splitext(meshfile)
	meshfile = meshfilebase * ".nas"
	mesh = import_NASTRAN("meshes/" * meshfile)
	fens = mesh["fens"]
	fes = mesh["fesets"][1]
    for i in 2:length(mesh["fesets"])
        fes = cat(fes, mesh["fesets"][i])
    end
    @show count(fes)

	r1 = selectelem(fens, fes, box = [-0.1 0.12 0 0.1 0.22 0.5], inflate = 1.0e-3)
	r2 = selectelem(fens, fes, box = [-0.1 0.12 0 0.1 0.0 0.2], inflate = 1.0e-3)
	r3 = setdiff(1:count(fes), vcat(r1, r2))
    @show count(fes), length(r1),length(r2), length(r3)

	fesets = [subset(fes, r1) subset(fes, r2) subset(fes, r3)]

	integrationrulestiff = TetRule(4)
	integrationrulemass = TetRule(4)
	
	return fens, fesets, integrationrulestiff, integrationrulemass
end

function lug_setup(E, nu, rho, fens, fesets, integrationrulestiff = nothing, integrationrulemass = nothing, massformulation = :lumped)
	@assert length(fesets) == 3

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
	fenids = selectnode(fens, plane = [1.0 0.0 0.0 0.22], inflate = 1.0e-5)
	setebc!(u, fenids, true, 1, 0.0)
	setebc!(u, fenids, true, 2, 0.0)
	setebc!(u, fenids, true, 3, 0.0)
	applyebc!(u)
	numberdofs!(u)
    @show count(fens), u.nfreedofs
	
	K = spzeros(u.nfreedofs, u.nfreedofs); M = spzeros(u.nfreedofs, u.nfreedofs)
	if (integrationrulestiff != nothing) && (integrationrulemass != nothing)
		MR = DeforModelRed3D
		material = MatDeforElastIso(MR, rho, E, nu, 0.0)

		for i = 1:length(fesets)
			femm = FEMMDeforLinear(MR, IntegDomain(fesets[i], integrationrulestiff), material)
			femm = associategeometry!(femm, geom)
			K  +=  stiffness(femm, geom, u)
		end
		K .= 0.5 * (K .+ transpose(K))

		massassembler = (massformulation == :lumped ? SysmatAssemblerSparseHRZLumpingSymm : SysmatAssemblerSparseSymm)
		for i = 1:length(fesets)
			femm = FEMMDeforLinear(MR, IntegDomain(fesets[i], integrationrulemass), material)
			femm = associategeometry!(femm, geom)
			M  += mass(femm, massassembler(), geom, u)
		end
		M .= 0.5 * (M .+ transpose(M))

		return u, K, M
	else
		return u
	end
end
