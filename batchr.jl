using Pkg; Pkg.activate("."); Pkg.instantiate()
using Revise

using FinEtoolsRapidHarmonicVA

@info "FinEtoolsRapidHarmonicVA ready"

include("examples/twisted_bar/define_sim.jl")

@info "Burn in "

# Burn in
for mesh_n in [6,]
    for nmodes in [50, ]
     
         sim2 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free", harmonic_method = "modal")
         solve(sim2, make_model)
         
         sim3 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced", harmonic_method = "modal")
         solve(sim3, make_model) 
         
         sim4 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "wyd_ritz", harmonic_method = "modal")
         solve(sim4, make_model) 
         
         sim5 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "reduced_conc", harmonic_method = "modal")
         solve(sim5, make_model) 

        #sim5 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced_enhanced", harmonic_method = "modal")
        #solve(sim5, make_model) 

    end
end

@info "Running simulations"

#Now for real
for mesh_n in [6, 12, 24]
    for nmodes in [50, ]

        sim2 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free", harmonic_method = "modal")
        solve(sim2, make_model)
        
        sim3 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced", harmonic_method = "modal")
        solve(sim3, make_model) 
        
        sim4 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "wyd_ritz", harmonic_method = "modal")
        solve(sim4, make_model) 
        
        sim5 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "reduced_conc", harmonic_method = "modal")
        solve(sim5, make_model) 

    end
end
