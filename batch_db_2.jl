using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

#using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/brake_disc/define_sim.jl")

#reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
reduction_methods = ["lanczos_ritz"]
#reduction_methods = ["free", "two_stage_free", "two_stage_wyd_ritz"]
reduction_methods = ["free", "wyd_ritz", "two_stage_free", "two_stage_wyd_ritz"]
reduction_methods = ["wyd_ritz", "two_stage_free", "two_stage_wyd_ritz"]
reduction_methods = ["wyd_ritz", "two_stage_free", "two_stage_wyd_ritz"]

@info "Burn in "

# Burn in

for nmodes in [50, ]

    for reduction_method in reduction_methods
        sim = define_sim(; nmodes = nmodes, reduction_method = reduction_method)
        #exportgraphics(sim_directory(), sim)
        solve(sim_directory(), sim, make_model)
    end
end


@info "Running simulations"

#Now for real

for nmodes in [50, 100,]

    for reduction_method in reduction_methods
        sim = define_sim(; nmodes = nmodes, reduction_method = reduction_method)
        #exportgraphics(sim_directory(), sim)
        solve(sim_directory(), sim, make_model)
    end
end
