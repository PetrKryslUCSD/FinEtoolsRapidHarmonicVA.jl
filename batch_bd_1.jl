using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

#using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/brake_disc/define_sim.jl")

reduction_methods = ["free", "two_stage_free_enh", "two_stage_free", "two_stage_wyd_ritz", "wyd_ritz"]
reduction_methods = ["free",]

@info "Burn in "

linsolve_method = "minres"
itmax = 5

for mesh_n in [1, ]
    for nmodes in [50, ]
        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
            for itmax in [5]
                # exportgraphics(sim_directory(), sim)
                solve(sim_directory(), sim, make_model)
            end
        end
    end
end


@info "Running simulations"

for mesh_n in [3, 4, 5]
    for nmodes in [100, 200, 400]
        for reduction_method in reduction_methods
            itmaxs = (reduction_method == "two_stage_free_enh") ? [5, 10, 20] : [0]
            for itmax in itmaxs
                sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
                    # exportgraphics(sim_directory(), sim)
                    solve(sim_directory(), sim, make_model)
            end
        end
    end
end
