using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

#reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
# reduction_methods = ["free",]
reduction_methods = ["free", "two_stage_free_residual", "two_stage_free", "two_stage_wyd_ritz"]
#reduction_methods = ["lanczos_ritz", "wyd_ritz"]
#reduction_methods = ["two_stage_free"]

@info "Burn in "

for mesh_n in [4, ]
   for nmodes in [25, ]

       simd = define_sim(;mesh_n=mesh_n, nmodes=0, reduction_method="none", harmonic_method="direct")
         # solve(sim_directory(), simd, make_model)

       sims = [simd]
       for reduction_method in reduction_methods
        sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
        solve(sim_directory(), sim, make_model)
        push!(sims, sim)
    end

end

end

@info "Running simulations"

#Now for real
for mesh_n in [4, 6, 8, 10, 12, 14, 16]
    for nmodes in [25, 50, 100, 200, 400]

        simd = define_sim(;mesh_n=mesh_n, nmodes=0, reduction_method="none", harmonic_method="direct")
        # solve(sim_directory(), simd, make_model)

        sims = [simd]
     
        for reduction_method in reduction_methods
            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
            solve(sim_directory(), sim, make_model)
            push!(sims, sim)
        end
        
    end
end
