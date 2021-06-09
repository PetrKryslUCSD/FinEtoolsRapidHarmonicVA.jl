using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

#using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/brake_disc/define_sim.jl")

reduction_methods = ["free", "two_stage_free", ]

@info "Burn in "

linsolve_method = "minres"
itmax = 5

# @info "Running simulations"

# for mesh_n in [3, ]
#     for nmodes in [100, ]
#         for reduction_method in reduction_methods
#             itmaxs = (reduction_method == "two_stage_free_enh") ? [20] : [0]
#             for itmax in itmaxs
#                 sim = define_sim(; namebase = "long", mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax )
#                 solve(sim_directory(), sim, make_model)
#             end
#         end
#     end
# end


reduction_methods = [ "two_stage_free_enh", ]

for mesh_n in [3, ]
    for nmodes in [100, ]
        for reduction_method in reduction_methods
            for itmax in [20]
                sim = define_sim(; namebase = "short", mesh_n = mesh_n, nmodes = nmodes, reduction_method = reduction_method, linsolve_method = linsolve_method, itmax = itmax, frequency_sweep = (1300, 3100, 800), resonance_list = 1:5)
                solve(sim_directory(), sim, make_model)
            end
        end
    end
end
