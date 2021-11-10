# Script to investigate the twisted bar example

include("./define_sim.jl")

reduction_methods = ["free", "wyd_ritz", "two_stage_free", "two_stage_wyd_ritz"]

@info "Burn in "

for mesh_n in [4, ]
    for nmodes in [50, ]
        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, )
            solve(sim_directory(), sim, make_model)
        end
    end
end

@info "Running simulations"

##Now for real

reduction_methods = ["free", "wyd_ritz", "two_stage_free", "two_stage_wyd_ritz"]

mesh_n = 4

for nmodes in [25, 50, ]
    for reduction_method in reduction_methods
        sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, )
        solve(sim_directory(), sim, make_model)
    end
end
