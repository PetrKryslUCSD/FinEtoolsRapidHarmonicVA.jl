import LinearAlgebra: eigen, qr, norm, mul!, dot, cholesky, Symmetric
using Statistics
using Arpack
using FinEtools
using FinEtoolsDeforLinear
import CoNCMOR: CoNCData, transfmatrix, LegendreBasis, SineCosineBasis 

function load_mesh(meshfile)
    output = import_MESH(meshfile)
    fens, fes = output["fens"], output["fesets"][1]
    return  fens, fes
end

function model_setup(sim)
    prop = retrieve_json(sim)

    meshfile = prop["mesh"]

    fens, fes = load_mesh(joinpath(prop["meshesdir"], meshfile))

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    applyebc!(u)
    numberdofs!(u)
    
    if isa(fes, FESetH8)
        integrationrulestiff = GaussRule(3,2)
        integrationrulemass = GaussRule(3,3)
    else
        integrationrulestiff = TetRule(1)
        integrationrulemass = TetRule(4)
    end

    E = prop["E"]
    nu = prop["nu"]
    rho = prop["rho"]

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, integrationrulestiff), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K .= 0.5 * (K .+ transpose(K))

    femm = FEMMDeforLinear(MR, IntegDomain(fes, integrationrulemass), material)
    femm = associategeometry!(femm, geom)
    M = mass(femm, SysmatAssemblerSparseHRZLumpingSymm(), geom, u)
    M .= 0.5 * (M .+ transpose(M))

    # Compute the parameters of Rayleigh damping. For the two selected
    # frequencies we have the relationship between the damping ratio and
    # the Rayleigh parameters
    a0(zeta1, zeta2, o1, o2) = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);# a0
    a1(zeta1, zeta2, o1, o2) = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);# a1

    zeta1, zeta2, omega1, omega2 = 0.001, 0.001, 2500, 7500
    C = a0(zeta1, zeta2, omega1, omega2) * M + a1(zeta1, zeta2, omega1, omega2) * K

    F = fill(0.0, size(K, 1))
    forcenl = selectnode(fens, nearestto=prop["force_location"])
    forcedir = prop["force_direction"]
    forcendof = u.dofnums[forcenl, forcedir]
    F[forcendof...] = 1000.0

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

    return model
end

function reduced_basis(sim)
    prop = retrieve_json(sim)
    if prop["reduction_method"] == "free_reduced"
        redu_free_vibration(sim)
    elseif prop["reduction_method"] == "free"
        full_free_vibration(sim)
    else
        @error "Unknown reduced-basis method"
    end
    true
end

function full_free_vibration(sim)
    @info "Free Vibration"
    prop = retrieve_json(sim)

    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

    results = Dict()
    if isfile(resultsfile)
        results = retrieve_json(resultsfile)
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = model_setup(sim)       
    end

    fens = model["fens"]
    @info "$(count(fens)) nodes"

    neigvs = prop["neigvs"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]
    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = eigs(Symmetric(K + mass_shift*M), Symmetric(M); nev=neigvs, which=:SM)
        eval .-=  mass_shift;
        fs = real(sqrt.(complex(eval)))/(2*pi)
    end
    println("Natural frequencies: $(round.(fs, digits=4)) [Hz]")

    rd = Dict()

    rd["frequencies"] = fs
    timing["Total"] = timing["Problem setup"] + timing["EV problem"]
    rd["timing"] = timing

    mkpath(prop["matricesdir"])
    rd["eigenvectors"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-evec", "h5")))
    store_matrix(rd["eigenvectors"]["file"], evec)
    rd["eigenvalues"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-eval", "h5")))
    store_matrix(rd["eigenvalues"]["file"], eval)

    results["free_vibration"] = rd
    store_json(resultsfile, results)

    true
end

function harmonic_vibration(sim)
    prop = retrieve_json(sim)
    if prop["harmonic_method"] == "modal"
        full_harmonic_vibration(sim)
    elseif prop["harmonic_method"] == "direct"
        full_harmonic_vibration_direct(sim)
    else
        @error "Unknown harmonic-vibration method"
    end
    true
end

function full_harmonic_vibration(sim)
    @info "Harmonic Vibration (modal)"
    prop = retrieve_json(sim)
    
    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

    if !isfile(resultsfile)
        @error "Need the results"
    end
    results = retrieve_json(resultsfile)

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = model_setup(sim)       
    end

    K = model["K"]
    M = model["M"]
    C = model["C"]
    F = model["F"]
    f = results["free_vibration"]["eigenvectors"]["file"]
    evecs = retrieve_matrix(f)
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
    f = results["free_vibration"]["eigenvalues"]["file"]
    evals = retrieve_matrix(f)
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

    mkpath(prop["matricesdir"])
    rd["frf"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-frf", "h5")))
    store_matrix(rd["frf"]["file"], frf)

    results["harmonic_vibration"] = rd
    store_json(resultsfile, results)

    true
end

function full_harmonic_vibration_direct(sim)
    @info "Harmonic Vibration (direct)"
    prop = retrieve_json(sim)
    
    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

    if !isfile(resultsfile)
        @error "Need the results"
    end
    results = retrieve_json(resultsfile)

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = model_setup(sim)       
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
    f = results["free_vibration"]["eigenvalues"]["file"]
    evals = retrieve_matrix(f)
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

    mkpath(prop["matricesdir"])
    rd["frf"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-frf", "h5")))
    store_matrix(rd["frf"]["file"], frf)

    results["harmonic_vibration"] = rd
    store_json(resultsfile, results)

    true
end


function reducedmodelparameters(V, N, E, nu, rho, fmax, alpha, smallestdimension = Inf)
    V, N, E, nu, rho, fmax, alpha
    c = sqrt(E / 2 / (1 + nu) / rho)
    lambda = c / fmax
    d = min(lambda, smallestdimension)
    Ncfloat = V / d^3
    Nchi = nextpow(2, Ncfloat)
    Nclo = prevpow(2, Ncfloat)
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
    nbf1max = nbf1max < 4 ? 4 : nbf1max
    nbf1max = nbf1max > 9 ? 9 : nbf1max
    return Nc, nbf1max
end

function redu_free_vibration(sim)
    @info "Free Vibration (Reduced)"
    prop = retrieve_json(sim)

    mkpath(prop["resultsdir"])
    rf = with_extension(sim * "-results", "json")
    resultsfile = joinpath(prop["resultsdir"], rf)

    results = Dict()
    if isfile(resultsfile)
        results = retrieve_json(resultsfile)
    end

    timing = Dict{String, FFlt}()

    timing["Problem setup"] = @elapsed begin
        model = model_setup(sim)       
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

    Nc, nbf1max = reducedmodelparameters(V, N, E, nu, rho, fmax, alpha)
    @info "Number of clusters $Nc, number of functions $nbf1max"

    timing["Partitioning"] = @elapsed begin
        partitioning = nodepartitioning(fens, Nc)
        mor = CoNCData(fens, partitioning)
    end

    u = model["u"]

    timing["Transformation matrix"] = @elapsed begin
        Phi = transfmatrix(mor, LegendreBasis, nbf1max, u);
    end

    neigvs = prop["neigvs"]
    mass_shift = prop["mass_shift"]
    K = model["K"]
    M = model["M"]

    transfm(m, evecs) = (evecs' * m * evecs)
    timing["Reduced matrices"] = @elapsed begin
        rK0 = transfm(K, Phi)
        rK0 .= 0.5 * (rK0 .+ transpose(rK0))
        rM0 = transfm(M, Phi)
        rM0 .= 0.5 * (rM0 .+ transpose(rM0))
    end
    @info "Transformation matrix dimensions $(size(Phi))"

    timing["EV problem"] = @elapsed begin
        eval, evec, nconv = eigs(rK0 + mass_shift*rM0, rM0; nev=neigvs, which=:SM)
        approxfs = @. real(sqrt(complex(eval - mass_shift)))/(2*pi);
        approxevec = Phi*real(evec)
    end

    println("Approximate natural frequencies: $(round.(approxfs, digits=4)) [Hz]")

    rd = Dict()

    rd["frequencies"] = approxfs
    timing["Total"] = timing["Problem setup"] + timing["Partitioning"] + timing["Transformation matrix"] + timing["Reduced matrices"] + timing["EV problem"]
    rd["timing"] = timing

    mkpath(prop["matricesdir"])
    rd["eigenvectors"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-evec", "h5")))
    store_matrix(rd["eigenvectors"]["file"], approxevec)
    rd["eigenvalues"] = Dict("file"=>joinpath(prop["matricesdir"], with_extension(sim * "-eval", "h5")))
    store_matrix(rd["eigenvalues"]["file"], eval)

    results["free_vibration"] = rd
    store_json(resultsfile, results)

    true
end

