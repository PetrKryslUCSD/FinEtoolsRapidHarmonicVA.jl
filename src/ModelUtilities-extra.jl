
function conc_reduced(cdir, sim, make_model)
    @info "CoNC (Reduced)"
    prop = retrieve_json(joinpath(cdir, sim))

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
    nbf1maxclamp = prop["nbf1maxclamp"]

    matricesdir = prop["matricesdir"]
    mkpath(joinpath(cdir, matricesdir))

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension, nbf1maxclamp)
    
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
    rd["basis"] = Dict("file"=>joinpath(matricesdir, with_extension(sim * "-Phi", "h5")))
    store_matrix(joinpath(cdir, rd["basis"]["file"]), Phi)

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
    nbf1maxclamp = prop["nbf1maxclamp"]

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension, nbf1maxclamp)
    
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
        eval, evec, nconv = _eigs(Kr + mass_shift*Mr, Mr, nmodes)
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
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

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

function two_stage_free_residual(cdir, sim, make_model)
    @info "Free Vibration (Reduced, Residual)"
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
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if "partitioning_method" in keys(model) && 
            model["partitioning_method"] == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            Nc = Int(round(N/1000))
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
            nbf1max = minimum(nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            V = integratefunction(femm, geom, (x) ->  1.0)
            alpha = prop["alpha"]
            smallestdimension = prop["smallestdimension"]
            nbf1maxclamp = prop["nbf1maxclamp"]
            Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension, nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
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
        
    # The following works, but it is expensive
    timing["Additional vectors"] = @elapsed begin
        @show frequencies = approxfs[1:2]
        @time  begin
        C = model["C"]
        Cr = transfm(C, Phi)
        Fr = transfv(F, Phi)
    end
    @time begin
        vs = []
        for f in frequencies
            omega = 2*pi*f;
            Ur = (-omega^2*Mr + (1im*omega)*Cr + Kr)\Fr;
            resu = diag((-omega^2*M + (1im*omega)*C + K), 0) .\ (F - (-omega^2*M + (1im*omega)*C + K)*(Phi*Ur))
            if (norm(imag.(resu)) > 1e-10)
                push!(vs, vec(imag.(resu)/norm(imag.(resu))))
            end
            if (norm(real.(resu)) > 1e-10)
                push!(vs, vec(real.(resu)/norm(real.(resu))))
            end
        end
    end

        # for k in axes(approxevec, 2)
        #     @show dot(approxevec[:, k], vs[1])
        # end
        
        # for k in 1:length(vs)
        #     approxevec = hcat(approxevec,  vs[k])
        # end
    end

    # norig = size(approxevec, 2) - length(vs)
    # timing["modifiedGS"] = @elapsed begin
    #     for i in 1:norig
    #         for v in norig+1:norig+length(vs)
    #             approxevec[:, v] .-= (dot(view(approxevec, :, v), view(approxevec, :, i)) / norm(view(approxevec, :, i))^2) * view(approxevec, :, i)
    #         end
    #     end
    #     for i in norig+1:norig+length(vs)-1
    #         for v in i+1:norig+length(vs)
    #             approxevec[:, v] .-= (dot(view(approxevec, :, v), view(approxevec, :, i)) / norm(view(approxevec, :, i))^2) * view(approxevec, :, i)
    #         end
    #     end
    #     invalid = falses(length(vs))
    #     ptr = norig+1
    #     for i in norig+1:norig+length(vs)
    #         @show i, norm(view(approxevec, :, i))
    #         if (norm(view(approxevec, :, i)) > 1e-14)
    #             (i != ptr) && (approxevec[:, ptr] .= approxevec[:, i])
    #             ptr += 1
    #         end
    #     end
    #     approxevec = approxevec[:, 1:ptr-1]
    # end

    timing["Eigenvector reconstruction"] = @elapsed begin
        approxevec = Phi*real(approxevec)
        for k in 1:length(vs)
            approxevec = hcat(approxevec,  vs[k])
        end
        # @show size(approxevec)
    end

    timing["orthogonalizeExtra"] = @elapsed begin
        norig = size(approxevec, 2) - length(vs)
        P = Matrix(approxevec[:, norig+1:end]'*approxevec[:, norig+1:end])
        eigenObj = eigen(P)
        selectVectors = findall(abs.(eigenObj.values) .> 1e-8)
        # @show selectVectors
        P = approxevec[:, norig+1:end] * eigenObj.vectors[:, selectVectors]
        # @show size(P)
        approxevec[:, norig+1:norig+size(P, 2)] .= P
        # @show size(approxevec)
        approxevec = approxevec[:, 1:norig+size(P, 2)]
    end
    @show size(approxevec)
    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = size(approxevec, 2)
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"] + timing["Additional vectors"] + timing["orthogonalizeExtra"]
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

function two_stage_free_enhanced(cdir, sim, make_model)
    @info "Free Vibration (Reduced, Enhanced)"
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
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if "partitioning_method" in keys(model) && 
            model["partitioning_method"] == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            Nc = Int(round(N/1000))
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
            nbf1max = minimum(nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            V = integratefunction(femm, geom, (x) ->  1.0)
            alpha = prop["alpha"]
            smallestdimension = prop["smallestdimension"]
            nbf1maxclamp = prop["nbf1maxclamp"]
            Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension, nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
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
        
    # The following works, but it is expensive
    timing["Additional vectors"] = @elapsed begin
        C = model["C"]
        Cr = transfm(C, Phi)
        Fr = transfv(F, Phi)
        na = 2
        vs = []
        for k in 1:na
            @show omega = 2*pi*approxfs[k];
            # Ur = (-omega^2*Mr + (1im*omega)*Cr + Kr)\Fr;
            Ur = (-omega^2*M + (1im*omega)*C + K)\F;
            push!(vs, imag.(Ur)/norm(imag.(Ur)))
            push!(vs, real.(Ur)/norm(real.(Ur)))
        end

        # for k in axes(approxevec, 2)
        #     @show dot(approxevec[:, k], vs[1])
        # end
        
        # for k in 1:length(vs)
        #     approxevec = hcat(approxevec,  vs[k])
        # end
    end

    # norig = size(approxevec, 2) - length(vs)
    # timing["modifiedGS"] = @elapsed begin
    #     for i in 1:norig
    #         for v in norig+1:norig+length(vs)
    #             approxevec[:, v] .-= (dot(view(approxevec, :, v), view(approxevec, :, i)) / norm(view(approxevec, :, i))^2) * view(approxevec, :, i)
    #         end
    #     end
    #     for i in norig+1:norig+length(vs)-1
    #         for v in i+1:norig+length(vs)
    #             approxevec[:, v] .-= (dot(view(approxevec, :, v), view(approxevec, :, i)) / norm(view(approxevec, :, i))^2) * view(approxevec, :, i)
    #         end
    #     end
    #     invalid = falses(length(vs))
    #     ptr = norig+1
    #     for i in norig+1:norig+length(vs)
    #         @show i, norm(view(approxevec, :, i))
    #         if (norm(view(approxevec, :, i)) > 1e-14)
    #             (i != ptr) && (approxevec[:, ptr] .= approxevec[:, i])
    #             ptr += 1
    #         end
    #     end
    #     approxevec = approxevec[:, 1:ptr-1]
    # end

    timing["Eigenvector reconstruction"] = @elapsed begin
        approxevec = Phi*real(approxevec)
        for k in 1:length(vs)
            approxevec = hcat(approxevec,  vs[k])
        end
    end
    @show size(approxevec)
    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

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

function conc_basis_only(cdir, sim, make_model)
    @info "conc basis only"
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
    timing["Partitioning"] = @elapsed begin
        partitioning = Int[]
        if "partitioning_method" in keys(model) && 
            model["partitioning_method"] == "metis"
            @info "Metis partitioning"
            C = connectionmatrix(femm, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            Nc = Int(round(N/1000))
            partitioning = Metis.partition(g, Nc; alg = :KWAY)
            nbf1max = minimum(nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
        else # Default: Recursive Inertial Bisection
            @info "RIB partitioning"
            V = integratefunction(femm, geom, (x) ->  1.0)
            alpha = prop["alpha"]
            smallestdimension = prop["smallestdimension"]
            nbf1maxclamp = prop["nbf1maxclamp"]
            Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension, nbf1maxclamp)
            @info "Number of clusters $Nc, number of functions $nbf1max"
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
    # timing["Reduced matrices"] = @elapsed begin
    #     Kr = transfm(K, Phi)
    #     #Kr .= 0.5 * (Kr .+ transpose(Kr))
    #     Mr = transfm(M, Phi)
    #     #Mr .= 0.5 * (Mr .+ transpose(Mr))
    # end
    @info "Transformation matrix dimensions $(size(Phi))"

    rd = Dict()

    rd["number_of_nodes"] = count(fens)
    rd["number_of_modes"] = nmodes
    rd["number_of_clusters"] = Nc
    rd["nbf1max"] = nbf1max

    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"]
    rd["timing"] = timing

    file = joinpath(matricesdir, with_extension(sim * "-Phi", "h5"))
    rd["basis"] = Dict("file"=>file)
    store_matrix(joinpath(cdir, rd["basis"]["file"]), Phi)
    file = joinpath(matricesdir, with_extension(sim * "-eval", "h5"))
   
    results["reduced_basis"] = rd
    store_json(joinpath(cdir, resultsfile), results)

    true
end
