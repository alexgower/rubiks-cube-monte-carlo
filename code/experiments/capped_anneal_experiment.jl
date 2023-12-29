# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using DelimitedFiles

include("../probes/capped_anneal.jl")




@inbounds @fastmath function capped_anneal_experiment(simulation_name::String, original_configuration::Vector{Matrix{Int64}}, energy_cap::Float64, swap_move_probability::Float64, T_1::Float64, T_0::Float64, N_T::Int64; verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false, infinite_temperature_steps::Int64=0)


    L = size(original_configuration[1],1)
    cube = RubiksCube(L)
    cube.configuration = deepcopy(original_configuration)

    energy_floor = energy(cube)
    println("Energy Floor: $energy_floor")

    # Cover everything in try/except clause
    try

        temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
        temperature_vector = vcat([100000.0 for i in 1:infinite_temperature_steps], temperature_vector)

        temperature_vector, E_by_temperature_step, correlation_function_by_temperature_step, energy_floor_configurations, energy_floor_configuration_correlation_function_values  =  capped_anneal!(cube, temperature_vector, energy_cap, energy_floor, swap_move_probability=swap_move_probability, verbose_annealing=verbose_annealing, verbose_metropolis_swap=verbose_metropolis_swap)

        println("Final Configuration:")
        println(cube.configuration)
        println("Final Energy")
        println(energy(cube))


            # Save Results ----------

        try

            # Write Energy Results to File ----------
            energy_results_filename = simulation_name * '_' * string(swap_move_probability) * "_energy_results.csv"

            touch(joinpath("results/capped_anneal_results",energy_results_filename))

            open(joinpath("results/capped_anneal_results",energy_results_filename), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, energy_cap=$energy_cap, energy_floor=$energy_floor \n")
                write(simulation_file, "Original Configuration = $original_configuration \n")
                write(simulation_file, "Temperature T, E(T) \n")
                
                for temperature_index in 1:N_T+1
                    write(simulation_file, "$(temperature_vector[temperature_index]), $(E_by_temperature_step[temperature_index]) \n")
                end

                write(simulation_file, "\n")
                write(simulation_file, "# Final Configuration: \n")
                write(simulation_file, "# $(cube.configuration) \n")
            end


            # Write Energy Floor Results to File ----------
            energy_floor_results_filename = simulation_name * '_' * string(swap_move_probability) * "_energy_floor_results.csv"

            touch(joinpath("results/capped_anneal_results",energy_floor_results_filename))

            open(joinpath("results/capped_anneal_results",energy_floor_results_filename), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, energy_cap=$energy_cap, energy_floor=$energy_floor \n")
                write(simulation_file, "Original Configuration = $original_configuration \n")
                write(simulation_file, "Energy Floor Configuration, Configuration Correlation Function Value \n")
                
                for configuration_index in eachindex(energy_floor_configurations)
                    write(simulation_file, "$(energy_floor_configurations[configuration_index]), $(energy_floor_configuration_correlation_function_values[configuration_index]) \n")
                end

                write(simulation_file, "\n")
                write(simulation_file, "# Final Configuration: \n")
                write(simulation_file, "# $(cube.configuration) \n")
            end
        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end
        

        
        try

            # Create energy plot ----------

            graph = plot([i for i in 1:length(temperature_vector)], E_by_temperature_step, xlabel="Time Step", ylabel="Energy", title="Rubik's Cube Capped Anneal, L=$L, N_T=$N_T", label="Energy")

            hline!(graph, [energy_floor], linestyle=:dash, color=:green, label="Energy Floor")
            hline!(graph, [energy_cap], linestyle=:dash, color=:blue, label="Energy Cap")

            savefig(graph, "results/capped_anneal_results/$simulation_name.png")

            # display(graph)

            # Create configuration correlation function plot ----------
            graph = plot([i for i in 1:length(temperature_vector)], correlation_function_by_temperature_step, xlabel="Time Step", ylabel="Configuration Correlation Function", title="Rubik's Cube Capped Anneal, L=$L, N_T=$N_T", label="Configuration Correlation Function")

            savefig(graph, "results/capped_anneal_results/$(simulation_name)_correlation_function.png")

            # display(graph)



        catch ex
            println("Cannot display or save results")
            showerror(stdout, ex)

        end

    catch ex
        println("General Error")
        showerror(stdout, ex)

        open("error_log.txt", "a") do error_file
            write(error_file, "Some Error: \n" * sprint(showerror(stdout, ex)) * "\n")
        end

    end

end


