using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

#reduction_methods = ["free", "two_stage_free", "wyd_ritz", "two_stage_wyd_ritz"]
# reduction_methods = ["free",]
#reduction_methods = ["free", "two_stage_free_enh", "two_stage_free", "two_stage_wyd_ritz"]
#reduction_methods = ["lanczos_ritz", "wyd_ritz"]
reduction_methods = ["two_stage_free_enh", "two_stage_free_enh"]

@info "Burn in "

for linsolve_method in ["minres"]
    for itmax in [5, 10, 20, 40]
        for mesh_n in [8, ]
            for nmodes in [50, ]
                namebase = "$(linsolve_method)_$(itmax)"
                for reduction_method in reduction_methods
                    sim = define_sim(; namebase = namebase, mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax)
                    solve(sim_directory(), sim, make_model)
                end
            end
        end
    end
end

#@info "Running simulations"

##Now for real
#for mesh_n in [4, 6, 8, 10, 12, 14, 16]
#    for nmodes in [25, 50, 100, 200, 400]

#        simd = define_sim(;mesh_n=mesh_n, nmodes=0, reduction_method="none", harmonic_method="direct")
#        # solve(sim_directory(), simd, make_model)

#        sims = [simd]
     
#        for reduction_method in reduction_methods
#            sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method)
#            solve(sim_directory(), sim, make_model)
#            push!(sims, sim)
#        end
        
#    end
#end
