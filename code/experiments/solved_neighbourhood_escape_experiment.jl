# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using DelimitedFiles

include("../probes/history_anneal.jl")





@inbounds @fastmath function solved_neighbourhood_escape_experiment(simulation_name::String, L::Int64, swap_move_probabilities::Vector{Float64}, T::Float64, N_steps::Int64, N_trials::Int64; verbose_metropolis_swap::Bool=false)


    # Cover everything in try/except clause so can print errors to file if running remotely
    try

        mean_normalized_energies_by_p_swap_and_step = zeros(length(swap_move_probabilities), N_steps+1)
        standard_deviation_normalized_energies_by_p_swap_and_step = zeros(length(swap_move_probabilities), N_steps+1)

        for (index,swap_move_probability) in pairs(swap_move_probabilities)
            simulation_name_to_use = "$(simulation_name)_T=$(T)_P_s=$(swap_move_probability)"

            all_trials_normalized_E_by_step = zeros(N_trials,N_steps+1)
            all_trials_normalized_E_by_step[:,1] .= -1.0 # Set initial energy to -1.0

            temperature_vector::Vector{Float64} = [T for t in 1:N_steps]

            # Iterate over trials
            @inbounds @simd for trial in 1:N_trials

                # Run Rubik's Cube Energy History Anneal ----------

                # Create Rubik's cube object and run annealing function on it
                cube = RubiksCube(L)

                _, this_trial_E_by_step =  history_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mix=false)

                all_trials_normalized_E_by_step[trial,2:N_steps+1] .= - this_trial_E_by_step ./ solved_configuration_energy(cube)

                

            end

            mean_normalized_energies_by_p_swap_and_step[index,:] .= mean(all_trials_normalized_E_by_step, dims=1)[:]
            standard_deviation_normalized_energies_by_p_swap_and_step[index,:] .= std(all_trials_normalized_E_by_step, dims=1)[:]

            # Save Results ----------
            try
                touch(joinpath("results/solved_neighbourhood_escape_results",simulation_name_to_use))

                open(joinpath("results/solved_neighbourhood_escape_results",simulation_name_to_use,), "w") do simulation_file
                    write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T=$T, N_steps=$N_steps \n")
                    write(simulation_file, "Step Index t, E(t) \n")
                    
                    for trial in N_trials
                        write(simulation_file, "Trial $trial \n")

                        for step in 1:N_steps
                            write(simulation_file, "$step, $(all_trials_normalized_E_by_step[trial,step]) \n")
                        end

                    end
                end

            catch ex
                
                println("Cannot save results to file")
                showerror(stdout, ex)

            end


            try
                # Create plot ----------

                # Create different plots for raw results        
                graph = plot([step for step in 1:N_steps+1], transpose(all_trials_normalized_E_by_step), xlabel="Steps Taken From Solved Configuration", ylabel="-Energy/Solved Energy", title="Solved Neighbourhood Escape, L=$L, P_swap=$(swap_move_probability), Trials=$N_trials, T=$T", titlefontsize=10, legend=false)    
                hline!(graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
                hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")

                # Save graph ----------
                savefig(graph, "results/solved_neighbourhood_escape_results/$(simulation_name_to_use).png")

    
            catch ex
    
                println("Cannot display or save plot")
                showerror(stdout, ex)
    
            end

        end

        try
        # Create combined plots for mean and standard deviation results----------
    
        mean_std_graph = plot([step for step in 1:N_steps+1], transpose(mean_normalized_energies_by_p_swap_and_step), yerr=transpose(standard_deviation_normalized_energies_by_p_swap_and_step), markerstrokecolor=:auto, xlabel="Steps Taken From Solved Configuration", ylabel="-Energy/Solved Energy", title="Solved Neighbourhood Escape, L=$L, P_swap=$(swap_move_probabilities[1]), Trials=$N_trials, T=$T", titlefontsize=10, labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))    
        hline!(mean_std_graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
        hline!(mean_std_graph, [-1.0], linestyle=:dash, color=:black, label="")

        mean_graph = plot([step for step in 1:N_steps+1], transpose(mean_normalized_energies_by_p_swap_and_step), xlabel="Steps Taken From Solved Configuration", ylabel="-Energy/Solved Energy", title="Solved Neighbourhood Escape, L=$L, Trials=$N_trials, T=$T", titlefontsize=10, labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))    
        hline!(mean_graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
        hline!(mean_graph, [-1.0], linestyle=:dash, color=:black, label="")

        # Save graph ----------
        savefig(mean_std_graph, "results/solved_neighbourhood_escape_results/$(simulation_name)_T=$(T)_mean_std.png")
        savefig(mean_graph, "results/solved_neighbourhood_escape_results/$(simulation_name)_T=$(T)__mean.png")


        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end    


    catch ex
        println("General Error")
        showerror(stdout, ex)

    end
end
