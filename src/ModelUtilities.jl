using Metis
import LinearAlgebra: eigen, qr, norm, mul!, dot, cholesky, Symmetric
using Statistics
using SparseArrays
using KrylovKit
using Krylov
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MatrixUtilityModule: export_sparse, import_sparse
using LinearAlgebra: eigen
import CoNCMOR: CoNCData, transfmatrix, LegendreBasis, SineCosineBasis 

# With Arpack 0.5, the explicit transform argument is needed to avoid
# unnecessary and expensive calculation.

_eigs(A, B, nmodes) = eigs(Symmetric(A), Symmetric(B); nev=nmodes, which=:SM, explicittransform=:none)

function load_mesh(meshfile)
    output = import_MESH(meshfile)
    fens, fes = output["fens"], output["fesets"][1]
    return  fens, fes
end

function ABAT!(C, A, B, AT, AB)
    M, N = size(A)
    @assert N == size(B, 1)
    @assert N == size(B, 2)
    mul!(AB, A, B)
    mul!(C, AB, AT)
    return C
end

function free(cdir, sim, make_model)
    @info "Free Vibration"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))
    
    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    F = model["F"]
    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = _eigs(K + mass_shift*M, M, nmodes)
        eval .-=  mass_shift;
        fs = real(sqrt.(complex(eval)))/(2*pi)
    end
    println("Natural frequencies: $(round.(fs, digits=4)) [Hz]")
    println("EV problem: $(round.(timing["EV problem"], digits=4)) [s]")

    ## Enhance with static solution instead of the last eigenvector
    #timing["Static vector"] = @elapsed begin
    #    evec = hcat(evec,  K\F)
    #end

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    rd["frequencies"] = fs
    timing["Total"] = timing["Problem setup"] + timing["EV problem"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), evec)
    file = joinpath(matricesdir, with_extension(sim * "-eval", "h5"))
    rd["eigenvalues"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["eigenvalues"]["file"]), eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, powertwo = true, smallestdimension = Inf, nbf1maxclamp = (3, 7))
    V, N, E, nu, rho, fmax, alpha
    c = sqrt(E / 2 / (1 + nu) / rho)
    lambda = c / fmax
    d = min(lambda, smallestdimension)
    Ncfloat = V / d^3
    if powertwo
        Nchi = nextpow(2, Ncfloat)
        Nclo =    Ncfloat > 1.0 ? prevpow(2, Ncfloat) : 2
    else
        Nchi = Int(round(Ncfloat))
        Nclo = Nchi
    end
    if Nchi - Ncfloat > Ncfloat - Nclo
        Nc = Nclo
    else
        Nc = Nchi
    end
    nbf1max = Int(floor((N / Nc)^(1.0/3))) 
    # The maximum number of 1d basis functions is reduced to assist with
    # arrangements of nodes that could prevent the columns of the
    # transformation matrix from being linearly independent.
    nbf1max = Int(floor(nbf1max / alpha))
    # We  clamp the number of basis functions $n_1$ to ensure both a
    # sufficient number of such functions to generate a reasonably rich
    # approximation, and to prevent a hugely expensive computation with too
    # many basis functions.
    nbf1max = nbf1max < nbf1maxclamp[1] ? nbf1maxclamp[1] : nbf1max
    nbf1max = nbf1max > nbf1maxclamp[2] ? nbf1maxclamp[2] : nbf1max
    return Nc, nbf1max
end

function two_stage_free(cdir, sim, make_model)
    @info "Two-stage Free Vibration"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]
    
    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    nbf1maxclamp = prop["nbf1maxclamp"]

    mor = nothing
    V = integratefunction(femm, geom, (x) ->  1.0)
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]
    nbf1maxclamp = prop["nbf1maxclamp"]
    partitioning_method = "partitioning_method" in keys(model) ?  model["partitioning_method"] : "rib"
    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, (partitioning_method == "rib"), smallestdimension, nbf1maxclamp)
    @info "Number of clusters $Nc, number of functions $nbf1max"
    
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if partitioning_method == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            partitioning = nodepartitioning(fens, Nc)
        end
        mor = CoNCData(fens, partitioning)
    end
    
    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end
    @info "Transformation matrix dimensions $(size(Phi))"
    @info "Sparsity of Phi: $(nnz(Phi)/prod(size(Phi)))"
    
    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    @info "Sparsity of K: $(nnz(K)/prod(size(K)))"
    @info "Sparsity of M: $(nnz(M)/prod(size(M)))"

    #timing["Reduced matrices"] = @elapsed begin
    #    PhiT = Matrix(Phi')
    #    AB = fill(zero(eltype(K)), size(Phi, 2), size(Phi, 1))
    #    Kr = fill(zero(eltype(K)), size(Phi, 2), size(Phi, 2))
    #    @time Kr = ABAT!(Kr, PhiT, K, Phi, AB)
    #    Mr = fill(zero(eltype(M)), size(Phi, 2), size(Phi, 2))
    #    @time Mr = ABAT!(Mr, PhiT, M, Phi, AB)
    #end
    transfm(m, t, tT) = (tT * m * t)
    transfv(v, t, tT) = (tT * v)
    timing["Reduced matrices"] = @elapsed begin
        PhiT = Phi'
        Kr = transfm(K, Phi, PhiT)
        #Kr .= 0.5 * (Kr .+ transpose(Kr))
        Mr = transfm(M, Phi, PhiT)
        #Mr .= 0.5 * (Mr .+ transpose(Mr))
    end

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = _eigs(Kr + mass_shift*Mr, Mr, nmodes)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
    end

    timing["Eigenvector reconstruction"] = @elapsed begin
        approxevec = Phi*real(evec)
    end

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), approxevec)
    file = joinpath(matricesdir, with_extension(sim * "-eval", "h5"))
    rd["eigenvalues"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["eigenvalues"]["file"]), eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function two_stage_free_enh(cdir, sim, make_model)
    @info "Free Vibration (Enhanced with partial solve)"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]
    V = integratefunction(femm, geom, (x) ->  1.0)

    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    nbf1maxclamp = prop["nbf1maxclamp"]

    mor = nothing
    V = integratefunction(femm, geom, (x) ->  1.0)
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]
    nbf1maxclamp = prop["nbf1maxclamp"]
    partitioning_method = "partitioning_method" in keys(model) ?  model["partitioning_method"] : "rib"
    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, (partitioning_method == "rib"), smallestdimension, nbf1maxclamp)
    @info "Number of clusters $Nc, number of functions $nbf1max"
    
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if partitioning_method == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            partitioning = nodepartitioning(fens, Nc)
        end
        mor = CoNCData(fens, partitioning)
    end

    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end

    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    F = model["F"]
    @info "Sparsity of K: $(nnz(K)/prod(size(K)))"
    @info "Sparsity of M: $(nnz(M)/prod(size(M)))"

    transfm(m, evecs) = (evecs' * m * evecs)
    transfv(v, evecs) = (evecs' * v)
    timing["Reduced matrices"] = @elapsed begin
        Kr = transfm(K, Phi)
        #Kr .= 0.5 * (Kr .+ transpose(Kr))
        Mr = transfm(M, Phi)
        #Mr .= 0.5 * (Mr .+ transpose(Mr))
    end
    @info "Transformation matrix dimensions $(size(Phi))"

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = _eigs(Kr + mass_shift*Mr, Mr, nmodes)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
    end
    approxevec = evec
    @info "EV problem: $(round(timing["EV problem"], digits=2))"

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")
    
    # Add additional vectors due to resonance residuals
    itmax = prop["itmax"]
    resonance_list = prop["resonance_list"]
    linsolve_method = prop["linsolve_method"]
    @info "Enhancement: $(linsolve_method) w $(itmax) iterations"
    timing["Additional vectors"] = @elapsed begin

        alg = nothing; met = nothing
        if linsolve_method == "minres"
            alg = MinresSolver(size(K, 1), size(K, 2), typeof(F))
            met = minres!
        elseif linsolve_method == "symmlq"
            alg = SymmlqSolver(size(K, 1), size(K, 2), typeof(F))
            met = symmlq!
        elseif linsolve_method == "diom"
            alg = DiomSolver(size(K, 1), size(K, 2), 20, typeof(F))
            met = diom!
        else
            @error "Unknown linear system solver"
        end

         # Reconstruct the approximate eigenvectors
        approxevec = Phi*real(approxevec)

        #dM = diag(M)
        #dK = diag(K)

        phi = fill(zero(eltype(approxevec)), size(Phi, 1))
        x0 = similar(phi)
        r0 = similar(phi)
        for (i, r) in enumerate(resonance_list)
            f = approxfs[r]
            omega = 2*pi*f;
            # Initial guess of the solution of the forced vibration at this
            # frequency. It is a good solution for the reduced system, but not
            # so good for the full system. The difference will be used as the
            # additional vector.
            phi .= view(approxevec, :, r)

            # Compute residual of the free vibration
            x0 .= phi
            r0 .= (-omega^2*(M*x0) + K*x0)

            # Without preconditioning
            (DU, stats) = met(alg, (-omega^2*M + K), r0; atol = 0.0, rtol = 0.0, itmax = itmax, verbose=1)
            
            # Normalize the vector to unit length
            DU /= norm(DU)

            # Replace the least significant modes: the eigenvectors for high
            # frequencies are replaced by the additional vectors.
            ni = size(approxevec, 2) - length(resonance_list) + i
            approxevec[:, ni] = DU

            # An alternative strategy is to append to the modes that came from
            # the eigenvalue analysis 
            #approxevec = hcat(approxevec, DU)
        end
   
    end
    @info "Additional vectors: $(round(timing["Additional vectors"], digits=2))"
     
    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = size(approxevec, 2)
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"] + timing["Additional vectors"] 
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), approxevec)
    file = joinpath(matricesdir, with_extension(sim * "-eval", "h5"))
    rd["eigenvalues"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["eigenvalues"]["file"]), eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function two_stage_free_resid(cdir, sim, make_model)
    @info "Free Vibration (Enhanced with residuals)"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]
    V = integratefunction(femm, geom, (x) ->  1.0)

    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    nbf1maxclamp = prop["nbf1maxclamp"]

    mor = nothing
    V = integratefunction(femm, geom, (x) ->  1.0)
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]
    nbf1maxclamp = prop["nbf1maxclamp"]
    partitioning_method = "partitioning_method" in keys(model) ?  model["partitioning_method"] : "rib"
    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, (partitioning_method == "rib"), smallestdimension, nbf1maxclamp)
    @info "Number of clusters $Nc, number of functions $nbf1max"
    
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if partitioning_method == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            partitioning = nodepartitioning(fens, Nc)
        end
        mor = CoNCData(fens, partitioning)
    end

    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end

    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    F = model["F"]
    @info "Sparsity of K: $(nnz(K)/prod(size(K)))"
    @info "Sparsity of M: $(nnz(M)/prod(size(M)))"

    transfm(m, evecs) = (evecs' * m * evecs)
    transfv(v, evecs) = (evecs' * v)
    timing["Reduced matrices"] = @elapsed begin
        Kr = transfm(K, Phi)
        #Kr .= 0.5 * (Kr .+ transpose(Kr))
        Mr = transfm(M, Phi)
        #Mr .= 0.5 * (Mr .+ transpose(Mr))
    end
    @info "Transformation matrix dimensions $(size(Phi))"

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = _eigs(Kr + mass_shift*Mr, Mr, nmodes)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
    end
    @info "EV problem: $(round(timing["EV problem"], digits=2))"

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")
    
    # Add additional vectors due to resonance residuals
    resonance_list = prop["resonance_list"]
    @info "Enhancement with residuals"
    timing["Additional vectors"] = @elapsed begin
        T = Phi*real(evec)
        norig = size(T, 2)
        Kr = transfm(K, T)
        Mr = transfm(M, T)
        C = model["C"]
        if C === nothing 
            C = K
            Cr = deepcopy(Kr)
        else
            Cr = transfm(C, evecs)
        end
        loss_factor = "loss_factor" in keys(model) ? model["loss_factor"] : 0.0
        Fr = transfv(F, T)
        approxevec = deepcopy(T)
        for (i, r) in enumerate(resonance_list)
            f = approxfs[r]
            omega = 2*pi*f;
            cmult = loss_factor == 0.0 ? omega : loss_factor
            Ur = (-omega^2*Mr + (1im*cmult)*Cr + Kr)\Fr;
            Kd = (-omega^2*M + (1im*cmult)*C + K)
            resid = F - Kd*(T*Ur)
            approxevec = hcat(approxevec, imag.(resid))
            approxevec = hcat(approxevec, real.(resid))
        end

        P = approxevec[:, norig+1:end]'*approxevec[:, norig+1:end]
        eigenObj = eigen(P)
        selectVectors = findall(abs.(eigenObj.values) .> 1e-8)
        P = approxevec[:, norig+1:end] * eigenObj.vectors[:, selectVectors]
        approxevec[:, norig+1:norig+size(P, 2)] .= P
        approxevec = approxevec[:, 1:norig+size(P, 2)]
    end
    @info "Additional vectors: $(round(timing["Additional vectors"], digits=2))"
     
    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = size(approxevec, 2)
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"] + timing["Additional vectors"] 
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), approxevec)
    file = joinpath(matricesdir, with_extension(sim * "-eval", "h5"))
    rd["eigenvalues"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["eigenvalues"]["file"]), eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function wyd_ritz(cdir, sim, make_model)
    @info "WYD Ritz"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))
    
    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]
    u = model["u"]

    nmodes = prop["nmodes"]
    K = model["K"]
    M = model["M"]
    F = model["F"]

    timing["Factorize stiffness"] = @elapsed begin
        # Factorize the stiffness matrix
        Kf = cholesky(K)
    end

    timing["Ritz-vector matrix"] = @elapsed begin
        # Refer for instance to Bayo, Wilson 1984, Joo et al 1989

        # Allocate the transformation matrix
        Phi = fill(0.0, size(K, 1), nmodes)
        p = fill(0.0, size(K, 1))
        q = fill(0.0, size(K, 1))

        # Solve for the static displacement
        p .= Kf \ F
        # Use it as the first basis vector
        #@show dot(p, M * p)
        #p ./= sqrt(dot(p, M * p))
        #@show dot(p, M * p)

        Phi[:, 1] .= p ./ sqrt(dot(p, M * p))
        
        # Generate the other basis vectors
        for i in 2:nmodes
            p .= Kf \ (M * view(Phi, :, i-1))
            q .= M * p
            for j in 1:i-1
                c = dot(q, view(Phi, :, j))
                p .-= c .* view(Phi, :, j)
            end
            Phi[:, i] .= p ./ sqrt(dot(p, M * p))
        end
        #temp =  M * Phi
        #for i in 1:nmodes, j in 1:nmodes
        #    @show i, j, dot(view(Phi, :, j), view(temp, :, i))
        #end
        #maxoff = 0.0
        #for i in 1:size(Phi, 2)
        #    for j in i:size(Phi, 2)
        #        m = dot(view(Phi, :, i), M * view(Phi, :, j))
        #        if i != j
        #            maxoff = max(maxoff, abs(m))
        #        end
        #        #@assert (i == j ? abs(1.0 - dot(view(Phi, :, i), M * view(Phi, :, j))) : abs(0.0 - dot(view(Phi, :, i), M * view(Phi, :, j)))) < 1.0e-6 
        #    end
        #end
        #@show maxoff
    end
    

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    timing["Total"] = timing["Problem setup"] + timing["Factorize stiffness"] + timing["Ritz-vector matrix"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), Phi)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function lanczos_ritz(cdir, sim, make_model)
    @info "Lanczos Ritz"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))
    
    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]
    u = model["u"]

    nmodes = prop["nmodes"]
    K = model["K"]
    M = model["M"]
    F = model["F"]

    timing["Factorize stiffness"] = @elapsed begin
        # Factorize the stiffness matrix
        Kf = cholesky(K)
    end

    timing["Ritz-vector matrix"] = @elapsed begin
        # Refer for instance to Nour-Omid, Clough 1984.
        # The implementation below follows Kim et al 2003.
        # The problem with the two-term recurrence is that it does
        # not manage to maintain mass-orthogonality of the generated vectors.

        # Allocate the transformation matrix
            Phi = fill(0.0, size(K, 1), nmodes)
            br = fill(0.0, size(K, 1))
            Ms = fill(0.0, size(K, 1))

            br .= Kf \ F
            mul!(Ms, M, br)
            beta = sqrt(dot(br, Ms))
            Phi[:, 1] .= br ./ beta
            
            for j in 1:nmodes-1
                qj = view(Phi, :, j)
                mul!(Ms, M, qj)
                br .= Kf \ Ms
                alpha = dot(br, Ms)
                br .-= alpha .* qj
                if j > 1
                    br .-= beta .* view(Phi, :, j-1)
                end
                mul!(Ms, M, br)
                beta = sqrt(dot(br, Ms))
                Phi[:, j+1] .= br ./ beta
                #for i in 1:j
                #    for k in i:j
                #        @show i, k, dot(view(Phi, :, i), M * view(Phi, :, k))
                #    end
                #end
            end
            
    end
    

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    timing["Total"] = timing["Problem setup"] + timing["Factorize stiffness"] + timing["Ritz-vector matrix"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), Phi)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function two_stage_wyd_ritz(cdir, sim, make_model)
    @info "Two-stage WYD Ritz"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))
    
    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    femm = model["femm"]
    geom = model["geom"]

    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    nbf1maxclamp = prop["nbf1maxclamp"]

    mor = nothing
    V = integratefunction(femm, geom, (x) ->  1.0)
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]
    nbf1maxclamp = prop["nbf1maxclamp"]
    partitioning_method = "partitioning_method" in keys(model) ?  model["partitioning_method"] : "rib"
    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, (partitioning_method == "rib"), smallestdimension, nbf1maxclamp)
    @info "Number of clusters $Nc, number of functions $nbf1max"
    
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if partitioning_method == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            partitioning = nodepartitioning(fens, Nc)
        end
        mor = CoNCData(fens, partitioning)
    end
    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end

    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    F = model["F"]

    transfm(m, evecs) = (evecs' * m * evecs)
    transfv(m, evecs) = (evecs' * m)
    timing["Reduced matrices"] = @elapsed begin
        Kr = transfm(K, Phi)
        Kr .= 0.5 * (Kr .+ transpose(Kr))
        Mr = transfm(M, Phi)
        Mr .= 0.5 * (Mr .+ transpose(Mr))
        Fr = transfv(F, Phi)
    end
    @info "Transformation matrix dimensions $(size(Phi))"


    timing["Factorize stiffness"] = @elapsed begin
        # Factorize the stiffness matrix
        Krf = cholesky(Kr)
    end

    timing["Ritz-vector matrix"] = @elapsed begin

        # Allocate the transformation matrix
        Phir = fill(0.0, size(Kr, 1), nmodes)
        p = fill(0.0, size(Kr, 1))
        q = fill(0.0, size(Kr, 1))

        # Solve for the static displacement
        p .= Krf \ Fr

        Phir[:, 1] .= p ./ sqrt(dot(p, Mr * p))
        
        # Generate the other basis vectors
        for i in 2:nmodes
            p .= Krf \ (Mr * view(Phir, :, i-1))
            q .= Mr * p
            for j in 1:i-1
                c = dot(q, view(Phir, :, j))
                p .-= c .* view(Phir, :, j)
            end
            Phir[:, i] .= p ./ sqrt(dot(p, Mr * p))
        end

        # Now reconstruct the full transformation matrix
        Phi = Phi*real(Phir)
        #temp =  M * Phi
        #for i in 1:nmodes, j in 1:nmodes
        #    @show i, j, dot(view(Phi, :, j), view(temp, :, i))
        #end
    end
    
    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["Factorize stiffness"] + timing["Ritz-vector matrix"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), Phi)
    
    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function reduced_basis(cdir, sim, make_model)
    prop = retrieve_json(joinpath(cdir, sim))
    if prop["reduction_method"] == "two_stage_free"
        two_stage_free(cdir, sim, make_model)
    elseif prop["reduction_method"] == "two_stage_free_enhanced"
        two_stage_free_enhanced(cdir, sim, make_model)
    elseif prop["reduction_method"] == "two_stage_free_resid"
        two_stage_free_resid(cdir, sim, make_model)
    elseif prop["reduction_method"] == "free"
        free(cdir, sim, make_model)
    elseif prop["reduction_method"] == "wyd_ritz"
        wyd_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "lanczos_ritz"
        lanczos_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "conc_reduced"
        conc_reduced(cdir, sim, make_model)
    elseif prop["reduction_method"] == "conc_basis_only"
        conc_basis_only(cdir, sim, make_model)
    elseif prop["reduction_method"] == "two_stage_wyd_ritz"
        two_stage_wyd_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "two_stage_free_enh"
        two_stage_free_enh(cdir, sim, make_model)       
    elseif prop["reduction_method"] == "none"
    else
        @error "Unknown reduced-basis method"
    end
    true
end

function harmonic_vibration_modal(cdir, sim, make_model)
    @info "Modal Harmonic Vibration"
    prop = retrieve_json(joinpath(cdir, sim))
    
    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    if !isfile(joinpath(cdir, resultsfile))
        @error "Need the results"
    end
    results = retrieve_json(joinpath(cdir, resultsfile))

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    K = model["K"]
    M = model["M"]
    C = model["C"]
    loss_factor = "loss_factor" in keys(model) ? model["loss_factor"] : 0.0
    F = model["F"]
    f = results["reduced_basis"]["basis"]["file"]
    evecs = retrieve_matrix(joinpath(cdir, f))
    @info "Number of modes: $(size(evecs, 2))"
        
    transfm(m, evecs) = (evecs' * m * evecs)
    transfv(v, evecs) = (evecs' * v)
    timing["Reduced matrices"] = @elapsed begin
        Mr = transfm(M, evecs)
        Kr = transfm(K, evecs)
        if C === nothing 
            Cr = deepcopy(Kr)
        else
            Cr = transfm(C, evecs)
        end
        Fr = transfv(F, evecs)
    end

    fens = model["fens"]
    u = model["u"]
    if "sensor_location" in keys( prop )
        sensornl = selectnode(fens, nearestto=prop["sensor_location"])
    else
        sensornl = model["sensor_node"]
    end
    sensordir = prop["sensor_direction"]
    sensorndof = u.dofnums[sensornl, sensordir]

    # Frequencies in the log space.
    frequencies = logspace(log10(prop["frequency_sweep"][1]), log10(prop["frequency_sweep"][2]), prop["frequency_sweep"][3])
    #f = results["reduced_basis"]["eigenvalues"]["file"]
    #evals = retrieve_matrix(f)
    frf = zeros(FCplxFlt, length(frequencies))
    
    timing["Frequency sweep"] = @elapsed begin
        U1 = zeros(FCplxFlt, u.nfreedofs)
        for k in 1:length(frequencies)
            omega = 2*pi*frequencies[k];
            # Solve the reduced equations. Is the damping viscous or structural?
            # If the loss factor is 0, viscous damping is assumed.
            cmult = loss_factor == 0.0 ? omega : loss_factor
            Ur = (-omega^2*Mr + 1im*cmult*Cr + Kr)\Fr;
            # Reconstruct the solution in the finite element space.
            U1 .= evecs * Ur;
            frf[k] = U1[sensorndof][1]
            # print(".")
        end
        # print("\n")
    end
    
    rd = Dict()

    rd["sweep_frequencies"] = frequencies
    timing["Total"] = timing["Problem setup"] + timing["Reduced matrices"] + timing["Frequency sweep"]
    rd["timing"] = timing

    file = joinpath(resultsdir, with_extension(sim * "-frf", "h5"))
    rd["frf"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["frf"]["file"]), frf)

    results["harmonic_vibration"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function harmonic_vibration_direct(cdir, sim, make_model)
    @info "Direct Harmonic Vibration"
    prop = retrieve_json(joinpath(cdir, sim))

    resultsdir = prop["resultsdir"]
    mkpath(joinpath(cdir, resultsdir))
    resultsfile = joinpath(resultsdir, with_extension(sim * "-results", "json"))
    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    results = Dict()
    if isfile(joinpath(cdir, resultsfile))
        results = retrieve_json(joinpath(cdir, resultsfile))
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = make_model(prop)       
    end

    K = model["K"]
    M = model["M"]
    C = model["C"]
    F = model["F"]
    
    fens = model["fens"]
    u = model["u"]
    sensornl = selectnode(fens, nearestto=prop["sensor_location"])
    sensordir = prop["sensor_direction"]
    sensorndof = u.dofnums[sensornl, sensordir]

    # Frequencies in the log space.
    frequencies = logspace(log10(prop["frequency_sweep"][1]), log10(prop["frequency_sweep"][2]), prop["frequency_sweep"][3])
    #f = results["reduced_basis"]["eigenvalues"]["file"]
    #evals = retrieve_matrix(f)
    frf = zeros(FCplxFlt, length(frequencies))
    
    timing["Frequency sweep"] = @elapsed begin
        U1 = zeros(FCplxFlt, u.nfreedofs)
        for k in 1:length(frequencies)
            omega = 2*pi*frequencies[k];
            # Solve the reduced equations.
            U1 .= (-omega^2*M + 1im*omega*C + K)\F;
            frf[k] = U1[sensorndof][1]
            print(".")
        end
        print("\n")
    end
    
    rd = Dict()

    rd["sweep_frequencies"] = frequencies
    timing["Total"] = timing["Problem setup"] + timing["Frequency sweep"]
    rd["timing"] = timing

    file = joinpath(resultsdir, with_extension(sim * "-frf", "h5"))
    rd["frf"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["frf"]["file"]), frf)

    results["harmonic_vibration"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function harmonic_vibration(cdir, sim, make_model)
    prop = retrieve_json(joinpath(cdir, sim))
    if prop["harmonic_method"] == "modal"
        harmonic_vibration_modal(cdir, sim, make_model)
    elseif prop["harmonic_method"] == "direct"
        harmonic_vibration_direct(cdir, sim, make_model)
    else
        @error "Unknown harmonic-vibration method"
    end
    true
end

function solve(cdir, sim, make_model)
    reduced_basis(cdir, sim, make_model)
    harmonic_vibration(cdir, sim, make_model)
    true
end