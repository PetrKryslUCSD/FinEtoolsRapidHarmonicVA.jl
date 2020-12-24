using FinEtoolsRapidHarmonicVA

include("sim_full.jl")
full_free_vibration(sim)
full_harmonic_vibration(sim)   
include("sim_redu.jl")
redu_free_vibration(sim)
#plot_mesh(sim)
full_harmonic_vibration(sim)   
#full_harmonic_vibration_direct(sim)    
plot_frf(["sim_full_nev_150", "sim_redu_nev_150" ])
plot_timing(["sim_full_nev_150", "sim_redu_nev_150" ])