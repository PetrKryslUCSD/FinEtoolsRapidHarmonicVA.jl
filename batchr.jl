using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

#reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
reduction_methods = ["lanczos_ritz", "wyd_ritz"]

@info "Burn in "

# Burn in
for mesh_n in [4, ]
    for nmodes in [25, ]
     
        for reduction_method in reduction_methods
           sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
           solve(sim_directory(), sim, make_model)
        end

    end
end

@info "Running simulations"

#Now for real
for mesh_n in [4, 6, ]
    for nmodes in [25, 50]

        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
            solve(sim_directory(), sim, make_model)
        end

    end
end
