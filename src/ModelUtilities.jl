import LinearAlgebra: eigen, qr, norm, mul!, dot, cholesky, Symmetric
using Statistics
using SparseArrays
using Arpack
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MatrixUtilityModule: export_sparse, import_sparse
import CoNCMOR: CoNCData, transfmatrix, LegendreBasis, SineCosineBasis 

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
        eval, evec, nconv = eigs(Symmetric(K + mass_shift*M), Symmetric(M); nev=nmodes, which=:SM)
        eval .-=  mass_shift;
        fs = real(sqrt.(complex(eval)))/(2*pi)
    end
    println("Natural frequencies: $(round.(fs, digits=4)) [Hz]")

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

function reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension = Inf)
    V, N, E, nu, rho, fmax, alpha
    c = sqrt(E / 2 / (1 + nu) / rho)
    lambda = c / fmax
    d = min(lambda, smallestdimension)
    Ncfloat = V / d^3
    Nchi = nextpow(2, Ncfloat)
    Nclo =    Ncfloat > 1.0 ? prevpow(2, Ncfloat) : 2
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
    nbf1max = nbf1max < 3 ? 3 : nbf1max
    nbf1max = nbf1max > 7 ? 7 : nbf1max
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
    V = integratefunction(femm, geom, (x) ->  1.0)

    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension)
    
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
    end

    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end
    @info "Transformation matrix dimensions $(size(Phi))"
    @info "Sparsity of Phi: $(nnz(Phi)/prod(size(Phi)))"
    #export_sparse("Phi.txt", Phi)

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
        eval, evec, nconv = eigs(Symmetric(Kr + mass_shift*Mr), Symmetric(Mr); nev=nmodes, which=:SM)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
        
    end

    timing["Eigenvector reconstruction"] = @elapsed begin
        approxevec = Phi*real(evec)
        #p = view(approxevec, :, 1)
        #approxevec[:, 1] .= p ./ sqrt(dot(p, M * p))
        #for i in 2:size(approxevec, 2)
        #    p = view(approxevec, :, i-1)
        #    q = M * p
        #    for j in i:size(approxevec, 2)
        #        c = dot(q, view(approxevec, :, j))
        #        approxevec[:, j] .-= c .* p
        #    end
        #    p = view(approxevec, :, i)
        #    approxevec[:, i] .= p ./ sqrt(dot(p, M * p))
        #end
        #for i in 1:size(approxevec, 2)
        #    for j in i:size(approxevec, 2)

        #        @assert (i == j ? abs(1.0 - dot(view(approxevec, :, i), M * view(approxevec, :, j))) : abs(0.0 - dot(view(approxevec, :, i), M * view(approxevec, :, j)))) < 1.0e-12 
        #    end
        #end
    end

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

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

function conc_reduced(sim, make_model)
    @info "CoNC (Reduced)"
    prop = retrieve_json(sim)

    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

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
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension)
    
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
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
    
    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] 
    rd["timing"] = timing

    mkpath(prop["matricesdir"])
    rd["basis"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-Phi", "h5")))
    store_matrix(rd["basis"]["file"], Phi)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function redu_free_vibration_alt(sim, make_model)
    @info "Free Vibration (Reduced)"
    prop = retrieve_json(sim)

    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

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
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension)
    
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
    end

    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end
    Phi = Matrix(Phi)
    @info "Transformation matrix dimensions $(size(Phi))"

    nmodes = prop["nmodes"]
    mass_shift = prop["mass_shift"]
     
    transfm(m, evecs) = (evecs' * m * evecs)
    timing["Reduced matrices"] = @elapsed begin
        # Direct assembly of the reduced matrices
        assembler = SysmatAssemblerReduced(Phi)
        Kr = stiffness(model["femm"], assembler, geom, u)
        Mr = mass(model["femm"], assembler, geom, u)
    end
    

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = eigs(Kr + mass_shift*Mr, Mr; nev=nmodes, which=:SM)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
        
    end

    timing["Mass orthogonalization"] = @elapsed begin
        approxevec = Phi*real(evec)
        #p = view(approxevec, :, 1)
        #approxevec[:, 1] .= p ./ sqrt(dot(p, M * p))
        #for i in 2:size(approxevec, 2)
        #    p = view(approxevec, :, i-1)
        #    q = M * p
        #    for j in i:size(approxevec, 2)
        #        c = dot(q, view(approxevec, :, j))
        #        approxevec[:, j] .-= c .* p
        #    end
        #    p = view(approxevec, :, i)
        #    approxevec[:, i] .= p ./ sqrt(dot(p, M * p))
        #end
        #for i in 1:size(approxevec, 2)
        #    for j in i:size(approxevec, 2)

        #        @assert (i == j ? abs(1.0 - dot(view(approxevec, :, i), M * view(approxevec, :, j))) : abs(0.0 - dot(view(approxevec, :, i), M * view(approxevec, :, j)))) < 1.0e-12 
        #    end
        #end
    end

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"]
    rd["timing"] = timing

    mkpath(prop["matricesdir"])
    rd["basis"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-Phi", "h5")))
    store_matrix(rd["basis"]["file"], approxevec)
    rd["eigenvalues"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-eval", "h5")))
    store_matrix(rd["eigenvalues"]["file"], eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end

function two_stage_free_enhanced(sim, make_model)
    @info "Free Vibration (Reduced, Enhanced)"
    prop = retrieve_json(sim)

    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

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
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension)
    
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
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
    transfv(v, evecs) = (evecs' * v)
    timing["Reduced matrices"] = @elapsed begin
        Kr = transfm(K, Phi)
        Kr .= 0.5 * (Kr .+ transpose(Kr))
        Mr = transfm(M, Phi)
        Mr .= 0.5 * (Mr .+ transpose(Mr))
    end
    @info "Transformation matrix dimensions $(size(Phi))"

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = eigs(Kr + mass_shift*Mr, Mr; nev=nmodes, which=:SM)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
        approxevec = Phi*real(evec)
    end
    
    # Enhance with static solution instead of the last eigenvector
    #timing["Additional vectors"] = @elapsed begin
    #    C = model["C"]
    #    transfm(m, evecs) = (evecs' * m * evecs)
    #    transfv(v, evecs) = (evecs' * v)
    #    Mr2 = transfm(M, approxevec)
    #    Kr2 = transfm(K, approxevec)
    #    Cr2 = transfm(C, approxevec)
    #    Fr2 = transfv(F, approxevec)
    #    Cd = spzeros(size(C, 1), size(C, 2))
    #    Kd = spzeros(size(C, 1), size(C, 2))
    #    for i in 1:size(C, 1)
    #        Cd[i, i] = C[i, i]
    #        Kd[i, i] = K[i, i]
    #    end
    #    na = 4
    #    vs = []
    #    for k in 1:na
    #        omega = 2*pi*approxfs[k];
            
    #        Ur2 = (-omega^2*Mr2 + (1im*omega)*Cr2 + Kr2)\Fr2;
    #        U1 = approxevec * Ur2;
    #        #Uk = (1/(-omega^2)) * (M \ (F - K * U1 -  (1im*omega)*(C * U1)))
    #        Uk = (-omega^2 * M + (1im*omega)*Cd + Kd) \ (F - K * U1 - (1im*omega)*(C * U1) + Kd * U1 + (1im*omega)*(Cd * U1))
    #        #Uk = (-omega^2*M + (1im*omega)*C + K) \ F
    #        v = real(Uk - U1)
    #        @show norm(v)
    #        v = v / norm(v)
    #        push!(vs, v)
    #        v = imag(Uk - U1)
    #        @show norm(v)
    #        v = v / norm(v)
    #        push!(vs, v)
    #    end
    #    for k in 1:length(vs)
    #        approxevec = hcat(approxevec,  vs[k])
    #    end
    #end
    
    #timing["Additional vectors"] = @elapsed begin
    #    C = model["C"]
    #    transfm(m, evecs) = (evecs' * m * evecs)
    #    transfv(v, evecs) = (evecs' * v)
    #    Mr2 = transfm(M, approxevec)
    #    Kr2 = transfm(K, approxevec)
    #    Cr2 = transfm(C, approxevec)
    #    Fr2 = transfv(F, approxevec)
    #    na = 4
    #    vs = []
    #    for k in 1:na
    #        omega = 2*pi*approxfs[k];
            
    #        Ur2 = (-omega^2*Mr2 + Kr2)\Fr2;
    #        U1 = approxevec * Ur2;
    #        #Uk = (1/(-omega^2)) * (M \ (F - K * U1 -  (1im*omega)*(C * U1)))
    #        Uk = (-omega^2 * M) \ (- K * U1)
    #        #Uk = (-omega^2*M + (1im*omega)*C + K) \ F
    #        v = real(Uk - U1)
    #        @show norm(v)
    #        v = v / norm(v)
    #        push!(vs, v)
    #    end
    #    for k in 1:length(vs)
    #        approxevec = hcat(approxevec,  vs[k])
    #    end
    #end

    # The following works, but it is expensive
    timing["Additional vectors"] = @elapsed begin
        C = model["C"]
        Cr = transfm(C, Phi)
        Fr = transfv(F, Phi)
        na = 2
        vs = []
        for k in 1:na
            @show omega = 2*pi*approxfs[k];
            Ur = (-omega^2*Mr + (1im*omega)*Cr + Kr)\Fr;
            U1 = Phi * Ur;
            #Uk = (1/(-omega^2)) * (M \ (F - K * U1 -  (1im*omega)*(C * U1)))
            Uk = (-omega^2*M + (1im*omega)*C + K) \ F
            v = real(Uk - U1)
            @show norm(v)
            v = v / norm(v)
            push!(vs, v)
            v = imag(Uk - U1)
            @show norm(v)
            v = v / norm(v)
            push!(vs, v)
        end
        for k in 1:length(vs)
            approxevec = hcat(approxevec,  vs[k])
        end
    end

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"]
    rd["timing"] = timing

    mkpath(prop["matricesdir"])
    rd["basis"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-Phi", "h5")))
    store_matrix(rd["basis"]["file"], approxevec)
    rd["eigenvalues"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-eval", "h5")))
    store_matrix(rd["eigenvalues"]["file"], eval)

    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

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
    F = model["F"]
    f = results["reduced_basis"]["basis"]["file"]
    evecs = retrieve_matrix(joinpath(cdir, f))
    @info "Number of modes: $(size(evecs, 2))"

    transfm(m, evecs) = (evecs' * m * evecs)
    transfv(v, evecs) = (evecs' * v)
    timing["Reduced matrices"] = @elapsed begin
        Mr = transfm(M, evecs)
        Kr = transfm(K, evecs)
        Cr = transfm(C, evecs)
        Fr = transfv(F, evecs)
    end

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
            Ur = (-omega^2*Mr + 1im*omega*Cr + Kr)\Fr;
            # Reconstruct the solution in the finite element space.
            U1 .= evecs * Ur;
            frf[k] = U1[sensorndof][1]
            print(".")
        end
        print("\n")
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
    V = integratefunction(femm, geom, (x) ->  1.0)

    N = count(fens)
    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]
    fmax = prop["fmax"]
    alpha = prop["alpha"]
    smallestdimension = prop["smallestdimension"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension)
    
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
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
    elseif prop["reduction_method"] == "free"
        free(cdir, sim, make_model)
    elseif prop["reduction_method"] == "wyd_ritz"
        wyd_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "lanczos_ritz"
        lanczos_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "conc_reduced"
        conc_reduced(cdir, sim, make_model)
    elseif prop["reduction_method"] == "two_stage_wyd_ritz"
        two_stage_wyd_ritz(cdir, sim, make_model)
    elseif prop["reduction_method"] == "none"
    else
        @error "Unknown reduced-basis method"
    end
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