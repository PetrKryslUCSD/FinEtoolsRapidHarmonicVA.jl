using FinEtoolsRapidHarmonicVA
include("define_sim.jl")

#sim1 = define_sim(; reduction_method = "free", harmonic_method = "modal")
#reduced_basis(sim1)
##plot_mesh(sim1)
#harmonic_vibration(sim1)   
##plot_frf([sim1, ])

#sim2 = define_sim(; reduction_method = "none", harmonic_method = "direct")
##reduced_basis(sim2)
###plot_mesh(sim)
##harmonic_vibration(sim2)

#sim3 = define_sim(; reduction_method = "free_reduced", harmonic_method = "modal")
#reduced_basis(sim3)
#harmonic_vibration(sim3)   

#sim4 = define_sim(; reduction_method = "wyd_ritz", harmonic_method = "modal")
#reduced_basis(sim4)
#harmonic_vibration(sim4)   


#sim5 = define_sim(; reduction_method = "free_reduced_enhanced", harmonic_method = "modal")
#reduced_basis(sim5)
#harmonic_vibration(sim5)   

using Plotting
using Plotting: plot_frf, plot_timing, plot_times_reduced_basis
using Plotting: plot_timing_reduced_basis, plot_frf_errors, plot_frf_amplitudes
using FinEtoolsRapidHarmonicVA

include("define_sim.jl")

using Plotting: correlation


#mesh_n = 8
#nmodes = 400

#sim1 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free", harmonic_method = "modal")
#sim2 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "wyd_ritz", harmonic_method = "modal")
#sim3 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced", harmonic_method = "modal")

#@show correlation(sim1, sim2)
#@show correlation(sim1, sim3)
#@show correlation(sim2, sim3)

for mesh_n in [4, ]
    for nmodes in [400, ]
     #
     sim0 = define_sim(; mesh_n = mesh_n, nmodes = 0, reduction_method = "none", harmonic_method = "direct")
     sim1 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free", harmonic_method = "modal")
     sim2 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "wyd_ritz", harmonic_method = "modal")
     sim3 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced", harmonic_method = "modal")
        
#        #sim5 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced_enhanced", harmonic_method = "modal")
#        #sim6 = define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "reduced_wyd_ritz", harmonic_method = "modal")
        #plot_frf([sim1, sim2, sim3], "frf-m$(mesh_n)-n$(nmodes).pdf")
        #plot_frf_errors([sim1, sim2, sim3,  sim4], "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        plot_frf_errors([sim0, sim1, sim2, sim3], "frf-errors-m$(mesh_n)-n$(nmodes).pdf")
        #plot_timing_reduced_basis([sim1, sim2, sim3, sim4], "timing-reduced_basis-m$(mesh_n)-n$(nmodes).pdf")
#        #plot_timing([sim1, sim3, sim4, sim5], "harmonic_vibration", "timing-harmonic_vibration-m$(mesh_n)-n$(nmodes).pdf")
    end
end


#mesh_n = 8
#list_of_sim_lists = []

#sim_list = []

#    for nmodes in [50, 100, 200, 400]
#        push!(sim_list, define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free", harmonic_method = "modal"))
#    end

#push!(list_of_sim_lists, sim_list)

#sim_list = []

#    for nmodes in [50, 100, 200, 400]
#        push!(sim_list, define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "wyd_ritz", harmonic_method = "modal"))
#    end

#push!(list_of_sim_lists, sim_list)

#sim_list = []

#    for nmodes in [50, 100, 200, 400]
#        push!(sim_list, define_sim(; mesh_n = mesh_n, nmodes = nmodes, reduction_method = "free_reduced", harmonic_method = "modal"))
#    end

#push!(list_of_sim_lists, sim_list)

#plot_times_reduced_basis(list_of_sim_lists)