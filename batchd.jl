using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

@info "Running simulations"

#Now for real
nmodes = 0
for mesh_n in [4, 6, 8]
        
        sim = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "none", harmonic_method = "direct")
        solve(sim_directory(), sim, make_model)
        
end
