using FinEtoolsRapidHarmonicVA

include("define_sim.jl")

sim1 = define_sim(; mesh_mm = 30, reduction_method = "free", harmonic_method = "modal")
reduced_basis(sim1)
#plot_mesh(sim)
harmonic_vibration(sim1)   

sim2 = define_sim(; mesh_mm = 30, reduction_method = "free", harmonic_method = "direct")
reduced_basis(sim2)
#plot_mesh(sim)
harmonic_vibration(sim2)

plot_frf([sim1, sim2])
plot_timing(["sim_full_mesh_15_nev_150", "sim_redu_mesh_15_nev_150" ], "free_vibration")
plot_timing(["sim_full_mesh_15_nev_150", "sim_redu_mesh_15_nev_150" ], "harmonic_vibration")