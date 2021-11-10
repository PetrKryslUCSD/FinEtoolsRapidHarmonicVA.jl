# Script to investigate the twisted bar example: direct formulation of the
# sweep

include("./define_sim.jl")

harmonic_method = "direct"
reduction_method = "none"

@info "Burn in "

for mesh_n in [4, ]
    sim = define_sim(; mesh_n = mesh_n, harmonic_method = harmonic_method, reduction_method = reduction_method)
    solve(sim_directory(), sim, make_model)
end

@info "Running simulations"

##Now for real

for mesh_n in [4, ]
    sim = define_sim(; mesh_n = mesh_n, harmonic_method = harmonic_method, reduction_method = reduction_method)
    solve(sim_directory(), sim, make_model)
end

