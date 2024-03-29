using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

#reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
# reduction_methods = ["free",]
reduction_methods = ["free", "two_stage_free_enh", "two_stage_free", "two_stage_wyd_ritz", "wyd_ritz"]
#reduction_methods = ["lanczos_ritz", "wyd_ritz"]
#reduction_methods = ["two_stage_free_enh", ]

@info "Burn in "

for linsolve_method in ["minres", "diom", "symmlq"]
    for itmax in [5]
        for mesh_n in [4, ]
            for nmodes in [50, ]
                for reduction_method in reduction_methods
                    sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
                    solve(sim_directory(), sim, make_model)
                end
            end
        end
    end
end

@info "Running simulations"

##Now for real

linsolve_method = "minres"
for itmax in [5, 10, 20, 40]
    for nmodes in [25, 50, 100, 200, 400]
        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
            solve(sim_directory(), sim, make_model)
        end
    end
end