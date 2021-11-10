# Script to exercise the example of the twisted bar with enhancement to the
# reduced basis.

include("./define_sim.jl")

reduction_methods = ["free", "two_stage_free_enh", "two_stage_free", "two_stage_wyd_ritz", "wyd_ritz"]

@info "Burn in "

for linsolve_method in ["minres", ]
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
mesh_n = 4
for itmax in [5, 10, ]
    for nmodes in [25, 50, ]
        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
            solve(sim_directory(), sim, make_model)
        end
    end
end