using Plots
using DelimitedFiles

include("../core/rubiks_cube.jl")
include("../core/monte_carlo.jl")
include("../probes/history_anneal.jl")

@inbounds @fastmath function time_experiment(cube::RubiksCube, simulation_name::String, L::Int64, swap_move_probabilities::Vector{Float64}, T::Float64, N_t::Int64; verbose_metropolis_swap::Bool=false, normalization::String="solved")

    # Cover everything in try/except clause so can print errors to file if running remotely
    try
        time_vector::Vector{Int64} = [t for t in 1:N_t]
        temperature_vector::Vector{Float64} = [T for t in 1:N_t]
        normalised_E_by_time = []
        infinite_temperature_normalised_E_by_time = []

        for (index,swap_move_probability) in pairs(swap_move_probabilities)
            simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)

            # Run Rubik's Cube Energy History Anneal ----------

            # Use provided Rubik's cube object and run annealing function on it

            temperature_vector, E_by_time =  history_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mix=false)

            println("Final Configuration:")
            println(cube.configuration)

            push!(normalised_E_by_time, -E_by_time ./ solved_configuration_energy(cube))
            push!(infinite_temperature_normalised_E_by_time, -E_by_time ./ infinite_temperature_energy(cube))


            # Save Results ----------
            try
                touch(joinpath("results/funnel_escape_results",simulation_name_to_use))

                open(joinpath("results/funnel_escape_results",simulation_name_to_use), "w") do simulation_file
                    write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T=$T, N_t=$N_t \n")
                    write(simulation_file, "Time Index t, E(t) \n")
                    
                    for time_index in time_vector
                        write(simulation_file, "$time_index, $(E_by_time[time_index]) \n")
                    end
                end

            catch ex
                
                println("Cannot save results to file")
                showerror(stdout, ex)

            end
        end
        

        
        try
            # Create plot ----------
            
            if normalization == "solved"
                graph = plot(time_vector, normalised_E_by_time, xlabel="Time Step", ylabel="-Energy/Solved Energy", title="Rubik's Cube History Anneal, L=$L, N_t=$N_t", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
                hline!(graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
                hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")
            else # normalization == "infinite_temperature" case
                graph = plot(time_vector, infinite_temperature_normalised_E_by_time, xlabel="Time Step", ylabel="-Energy/Infinite Temperature Energy", title="Rubik's Cube History Anneal, L=$L, N_t=$N_t", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
                hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")
                hline!(graph, [-6.0], linestyle=:dash, color=:black, label="")
            end


            # Add other data to graph for comparison if exists ----------
            if isfile(joinpath("results/funnel_escape_results","other_data.csv"))
                data_matrix = readdlm(joinpath("results/funnel_escape_results","other_data.csv"), ',', Float64, '\n', skipstart=0)

                other_temperature_vector = copy(data_matrix[:,1])
                other_data = data_matrix[:,2]

                other_data = - other_data ./ solved_configuration_energy(RubiksCube(L))
            
                plot!(graph, other_temperature_vector, other_data, label="Other Results", seriestype=:scatter, color="blue", ms=2, ma=0.5)
            end

            # Save graph ----------
            savefig(graph, "results/funnel_escape_results/$simulation_name.png")

        catch ex

            println("Cannot display or save results")
            showerror(stdout, ex)

        end

    catch ex
        println("General Error")
        showerror(stdout, ex)

    end
end









function thin_funnel_experiment(N_t::Int64; normalization::String="solved")
    E = 2.5
    T = 2.1
    p_swap = 0.0

    # TODO redo E=2 for completion?

    # E_3_starting_configurations = [ [[3 4 4 4 4 6 6 5 2; 6 4 4 4 4 4 6 4 5; 4 2 4 4 2 2 5 5 5; 6 5 4 4 1 3 5 4 3; 4 5 4 1 1 1 5 4 4; 4 5 5 4 4 1 3 4 4; 4 5 3 5 5 5 3 3 2; 2 5 3 1 1 1 3 2 2; 5 5 2 1 1 1 2 2 4], [4 4 6 6 5 4 6 6 1; 4 2 4 6 4 6 6 6 4; 4 2 2 6 4 6 6 6 3; 4 4 4 4 4 2 2 3 3; 2 1 1 4 2 2 3 3 3; 2 6 6 5 2 2 2 1 1; 6 6 6 6 4 1 1 4 1; 6 3 1 5 5 1 1 1 1; 6 3 1 3 3 1 1 1 1], [2 5 4 5 5 5 5 2 3; 4 1 1 2 6 6 4 3 3; 4 1 1 1 6 2 4 3 1; 1 1 1 1 1 2 6 3 1; 1 1 1 5 3 2 2 3 3; 5 5 1 1 3 3 3 3 3; 5 3 1 1 3 3 3 3 3; 5 5 5 6 6 5 5 3 3; 5 5 1 2 1 5 2 3 2], [1 3 3 6 4 5 3 5 6; 1 6 6 6 6 5 5 5 5; 5 6 6 6 6 5 5 5 5; 3 3 3 6 6 5 4 4 5; 6 3 3 3 4 4 4 4 6; 3 3 3 3 3 4 4 4 6; 3 3 3 3 3 3 6 1 6; 3 3 3 3 3 3 6 6 6; 3 3 3 3 3 3 1 6 6], [3 1 1 1 5 5 4 1 4; 2 2 1 1 5 5 4 1 1; 2 2 2 2 5 1 4 4 3; 2 2 4 6 3 1 1 1 2; 2 2 2 5 5 5 1 1 5; 2 2 2 5 5 3 6 1 5; 5 2 2 1 1 3 5 5 5; 6 6 6 2 2 2 2 1 3; 5 4 6 2 2 2 2 4 4], [5 2 2 4 2 6 5 6 6; 4 4 4 4 2 3 5 5 6; 3 1 1 4 6 5 5 1 6; 4 2 5 6 6 6 5 6 6; 6 5 5 6 6 6 6 6 6; 6 5 2 2 2 5 6 6 2; 6 2 2 2 2 4 4 4 2; 1 2 2 2 2 2 2 4 2; 1 1 1 1 1 4 4 2 2]],
    #                             [[1 1 1 1 1 1 2 1 1; 6 1 1 1 6 6 1 1 1; 4 2 2 1 6 6 2 1 1; 5 5 2 4 1 1 1 1 1; 5 5 5 5 1 1 1 1 2; 5 5 5 1 1 1 1 2 4; 5 5 1 3 3 1 5 2 2; 5 3 3 3 3 3 3 2 2; 6 1 3 3 3 2 2 2 6], [5 5 5 2 4 4 4 3 4; 5 5 5 2 4 4 4 4 4; 5 5 2 2 4 4 4 4 4; 2 1 1 1 2 2 1 4 4; 3 1 3 3 2 6 6 4 3; 3 3 3 3 5 6 6 6 6; 3 3 5 5 5 4 6 6 6; 6 6 6 4 4 4 6 6 6; 6 6 6 4 4 4 6 4 6], [4 6 1 3 3 3 1 5 5; 6 1 4 4 3 3 3 5 3; 5 5 5 4 3 3 3 2 2; 5 5 5 2 3 3 3 4 5; 2 5 5 5 3 3 3 5 5; 2 6 6 3 3 3 3 1 5; 2 2 6 2 4 3 4 1 1; 2 2 6 3 3 2 4 2 3; 2 4 4 6 6 6 6 6 3], [3 5 6 5 5 3 3 3 3; 3 6 6 6 1 1 4 3 3; 1 6 6 6 1 1 1 3 4; 1 6 6 6 2 6 5 5 4; 1 1 4 4 4 4 4 4 4; 6 1 4 4 4 5 3 4 4; 1 1 4 4 1 5 4 4 2; 1 1 4 2 2 6 2 4 4; 1 1 4 2 2 6 4 2 2], [3 1 3 1 4 6 6 2 2; 1 3 3 3 6 1 1 4 2; 1 3 3 2 2 1 1 4 3; 1 4 2 2 2 2 2 2 2; 1 3 2 2 5 4 2 2 2; 3 3 2 5 1 4 4 2 2; 2 1 3 2 1 4 2 5 3; 3 2 2 2 2 5 5 5 2; 5 2 2 2 6 6 6 4 4], [4 4 4 1 6 6 5 5 5; 4 4 6 6 6 5 5 5 5; 5 2 6 6 2 5 5 2 5; 5 6 6 6 6 5 5 2 4; 6 6 6 6 6 5 5 5 1; 3 1 3 4 6 5 5 5 1; 5 6 3 4 6 6 1 1 6; 3 6 3 3 2 5 5 3 5; 2 4 3 3 5 5 3 6 1]],
    #                             [[5 5 3 3 1 3 3 3 3; 5 5 3 3 3 3 3 1 1; 4 3 3 3 3 3 3 5 5; 1 1 1 1 1 3 3 5 5; 2 1 5 5 1 3 3 2 2; 2 2 5 5 1 3 3 3 3; 2 5 5 3 3 3 3 1 3; 1 4 4 4 3 3 6 1 1; 1 4 4 2 2 1 4 1 1], [4 3 1 3 2 6 6 2 2; 4 3 3 3 2 6 4 4 2; 5 5 4 4 4 6 2 4 1; 5 5 4 4 5 3 3 2 2; 5 5 5 4 2 2 2 2 4; 5 5 5 1 4 4 1 5 1; 2 5 5 5 1 4 1 2 2; 2 2 1 2 1 6 1 2 2; 2 2 2 2 5 6 6 4 4], [2 2 1 2 3 4 1 4 2; 2 2 2 2 5 3 1 1 2; 2 2 2 2 3 3 6 1 2; 2 2 2 2 2 2 2 3 4; 4 2 2 3 3 2 6 6 4; 4 6 6 3 6 6 6 6 2; 5 6 6 2 2 5 3 3 5; 5 6 6 4 4 4 3 3 3; 5 6 6 4 6 3 3 3 3], [6 5 6 6 1 6 6 5 5; 5 5 5 5 5 5 5 5 5; 5 5 5 1 5 5 5 2 5; 5 1 5 1 1 1 4 1 5; 6 4 1 1 4 4 4 4 3; 4 4 1 4 4 4 4 4 4; 4 4 1 4 6 4 2 6 4; 4 5 5 3 6 4 4 3 6; 5 5 5 3 6 6 6 6 6], [3 3 3 1 3 3 3 3 3; 3 3 1 1 3 2 2 2 4; 3 3 1 1 2 2 2 2 1; 5 1 1 2 2 2 2 4 2; 5 5 5 5 5 3 1 4 3; 5 5 5 5 5 5 5 6 6; 5 4 4 2 4 2 4 4 1; 1 4 4 4 1 2 2 4 4; 1 1 1 1 1 3 4 3 4], [6 6 6 6 6 5 2 6 6; 6 6 6 6 6 5 2 6 6; 4 6 6 6 6 6 1 1 2; 4 6 6 6 6 5 1 1 1; 5 6 6 6 6 6 1 1 1; 4 1 6 6 3 6 1 2 6; 4 3 4 4 4 6 6 6 6; 6 1 1 1 3 6 6 6 1; 4 1 1 1 4 1 3 4 1]],
    #                             [[1 2 5 5 4 4 4 2 2; 1 2 2 5 5 3 4 2 4; 2 2 2 2 1 5 4 4 4; 3 5 2 2 1 5 4 4 1; 6 2 2 2 1 4 4 4 1; 5 2 2 5 6 6 5 3 3; 6 4 1 1 3 3 3 6 3; 6 4 1 1 5 6 6 6 3; 6 4 3 1 5 4 6 6 3], [5 1 3 3 3 1 1 5 5; 6 4 4 4 2 2 3 2 1; 6 6 6 5 2 6 2 1 1; 3 1 5 6 2 2 2 2 3; 1 1 1 5 2 2 1 2 3; 1 1 1 1 1 1 1 4 4; 1 1 1 1 1 1 1 4 4; 1 1 1 1 1 1 1 3 3; 2 1 1 1 1 1 1 3 3], [4 6 2 2 2 2 2 2 2; 2 2 2 6 3 3 5 5 2; 2 2 2 6 3 3 3 3 4; 2 2 3 3 3 3 3 3 4; 3 2 5 3 3 3 3 3 3; 3 3 5 3 3 3 3 3 3; 4 4 6 4 4 5 5 3 3; 4 4 4 4 4 5 2 3 3; 4 4 4 4 4 4 2 4 3], [4 5 4 2 4 6 3 3 4; 4 6 3 4 4 2 2 4 2; 4 6 1 6 5 4 2 2 2; 1 6 2 2 4 4 4 2 2; 1 3 5 4 4 1 2 6 2; 3 3 3 4 5 5 5 5 2; 3 3 4 4 5 5 5 3 2; 3 5 2 2 1 5 5 5 2; 3 4 6 5 5 5 5 5 6], [1 3 3 4 2 6 5 5 5; 1 3 3 3 3 1 1 1 5; 5 5 3 3 3 3 3 1 3; 5 4 2 2 1 1 1 1 2; 5 5 4 2 5 5 2 1 4; 6 5 4 4 4 4 2 5 4; 6 5 4 4 4 4 4 4 1; 6 1 5 4 4 4 6 3 1; 6 5 5 2 6 6 6 3 1], [6 2 5 6 6 6 6 6 2; 6 6 6 2 6 6 5 6 4; 6 3 6 2 6 6 6 6 2; 5 6 1 5 6 6 6 6 6; 5 5 6 5 6 6 6 6 2; 5 6 1 1 6 6 6 5 5; 5 5 5 6 6 6 5 1 1; 5 5 6 6 6 1 5 1 6; 5 5 5 6 6 1 1 1 1]],
    #                             [[5 2 2 4 3 4 4 4 4; 5 3 3 2 2 4 4 4 4; 1 3 3 3 6 4 4 4 3; 1 3 4 1 1 1 1 4 3; 1 1 1 1 1 1 1 1 4; 1 1 1 2 2 3 5 5 4; 1 1 1 1 1 5 6 6 6; 1 1 1 1 1 1 6 6 6; 1 1 1 1 1 1 6 6 6], [4 5 5 3 3 3 3 3 3; 4 5 4 3 3 3 3 3 3; 5 4 5 3 3 6 6 3 3; 5 4 2 6 1 6 3 3 3; 4 2 1 5 2 2 3 3 3; 4 6 6 3 2 2 1 3 3; 3 3 3 3 2 2 4 4 3; 3 2 2 2 2 1 1 5 3; 1 1 1 5 2 6 6 3 3], [2 2 6 3 2 2 5 5 6; 6 1 4 5 6 5 5 2 4; 4 5 5 6 6 2 2 2 4; 1 1 2 4 6 2 2 5 2; 1 2 2 2 3 3 3 3 2; 1 2 6 1 3 3 5 3 2; 1 1 1 1 5 5 5 1 2; 3 1 1 1 4 2 2 2 1; 3 2 4 4 4 2 4 2 2], [4 6 4 4 4 4 2 2 2; 4 4 4 3 4 2 3 1 1; 1 4 4 4 4 4 3 3 3; 3 4 4 4 4 1 5 4 2; 3 4 4 4 4 4 2 5 2; 4 4 4 4 4 2 2 1 2; 2 2 2 5 4 2 2 2 2; 6 2 5 5 5 4 5 4 4; 6 6 3 2 5 5 4 4 4], [6 6 2 2 5 6 6 2 2; 6 6 2 2 6 6 6 4 2; 6 6 2 2 2 1 4 6 6; 6 2 3 5 5 5 3 6 6; 6 4 5 5 5 5 6 5 5; 5 5 6 6 3 5 3 6 5; 1 5 6 6 5 5 3 2 5; 1 3 3 6 6 6 6 6 5; 3 3 3 3 6 6 5 5 5], [5 5 5 5 5 5 5 5 5; 1 5 5 5 5 5 5 5 5; 2 1 1 4 5 5 5 5 5; 6 6 6 4 6 5 1 3 6; 6 6 6 6 6 3 3 3 6; 5 6 4 6 6 3 3 2 6; 2 6 6 6 4 1 1 2 6; 4 6 6 4 1 1 1 3 3; 1 2 4 1 1 1 1 1 1]],
    #                             [[5 6 6 4 2 2 2 2 3; 6 5 5 3 2 2 4 4 4; 6 4 3 3 4 4 4 4 4; 6 6 5 5 1 2 4 5 1; 5 3 1 2 1 2 4 2 6; 5 2 2 2 2 2 6 6 6; 1 2 2 6 5 1 6 5 5; 2 2 2 3 3 4 5 6 6; 2 6 2 4 4 6 4 2 4], [5 5 2 3 3 3 3 3 6; 5 5 2 5 5 5 3 3 3; 5 2 5 5 5 3 3 3 2; 2 2 5 1 5 3 3 2 2; 2 2 2 2 2 4 2 1 3; 5 1 6 6 5 4 2 6 4; 5 6 6 2 2 2 2 2 4; 1 4 6 5 5 5 5 6 3; 1 1 1 5 5 5 5 6 4], [5 3 3 3 2 1 1 1 5; 5 3 3 6 6 1 1 1 3; 3 3 3 6 6 1 1 1 3; 3 3 3 3 1 1 1 1 1; 3 3 3 3 3 3 1 1 1; 3 3 3 3 3 6 1 1 1; 1 1 5 2 3 1 1 1 1; 1 1 1 3 3 6 6 1 1; 3 1 6 5 5 6 6 5 1], [3 1 1 1 1 1 3 3 6; 4 1 1 1 1 1 3 3 6; 2 5 1 1 1 1 4 3 2; 2 4 1 1 1 3 4 3 2; 1 4 4 4 4 4 4 5 2; 4 4 4 4 4 4 4 2 2; 3 4 4 4 1 4 4 3 4; 5 4 4 4 1 4 1 6 4; 4 4 4 6 6 4 4 4 3], [4 3 3 3 1 1 1 2 2; 4 4 4 1 4 1 1 2 2; 4 4 3 2 2 3 1 4 4; 4 4 4 4 6 5 2 4 4; 4 4 5 5 5 3 3 2 4; 5 3 5 5 5 2 2 2 4; 5 3 5 5 5 3 6 6 6; 1 3 5 3 6 6 6 6 5; 1 2 2 3 6 6 6 6 6], [6 4 6 6 6 5 3 2 2; 5 5 6 6 6 4 6 2 2; 5 5 6 6 6 3 2 2 2; 5 5 6 6 6 1 5 2 2; 4 4 6 6 6 1 3 5 5; 2 2 6 6 6 5 5 5 6; 1 2 2 6 6 5 5 6 6; 4 5 2 6 6 5 5 2 6; 1 3 5 3 3 1 5 5 2]],
    #                             [[5 6 6 4 4 3 3 3 3; 5 5 6 4 4 4 6 3 3; 1 1 3 6 4 2 2 2 3; 1 1 3 5 5 4 4 5 5; 1 1 3 1 1 2 5 5 5; 2 3 3 1 1 3 3 3 4; 2 4 4 1 1 4 6 6 4; 5 4 4 5 6 6 6 6 6; 5 4 3 1 5 5 6 4 4], [6 6 6 6 6 6 6 5 1; 6 6 6 6 6 6 6 4 4; 4 3 6 2 1 1 4 4 4; 4 1 6 1 1 4 4 4 4; 4 2 6 6 2 4 4 4 4; 2 1 2 2 2 4 4 4 4; 2 5 3 5 5 1 4 4 4; 4 5 1 5 5 4 4 4 4; 1 1 4 3 5 5 5 4 4], [6 1 1 1 2 2 5 5 5; 4 1 1 1 1 2 5 5 5; 5 1 1 1 2 2 1 3 4; 5 5 3 1 2 2 1 2 4; 3 3 3 3 3 5 1 1 1; 3 3 3 3 3 3 4 1 1; 3 3 4 4 3 3 2 3 1; 3 3 4 3 3 3 3 3 4; 3 3 3 3 3 3 1 3 4], [3 3 4 2 2 4 1 6 6; 3 4 5 6 2 4 4 6 6; 2 2 1 5 5 2 3 3 6; 2 2 1 2 4 3 3 3 3; 2 2 4 2 4 3 1 3 3; 2 2 2 2 4 5 6 6 1; 2 2 2 2 2 2 2 1 1; 2 2 2 4 4 2 1 1 1; 2 2 2 6 6 1 1 1 4], [2 2 5 6 6 5 5 1 2; 2 2 5 2 6 2 2 5 2; 2 2 5 5 6 5 1 4 2; 5 1 1 5 5 5 6 6 6; 4 1 2 5 5 4 4 6 6; 1 1 1 1 1 4 4 3 3; 3 1 5 5 2 4 3 3 3; 1 1 3 2 2 4 5 2 2; 1 1 6 6 1 6 5 2 2], [6 6 3 3 1 5 5 5 5; 6 3 2 5 5 5 5 2 5; 6 5 5 5 5 6 6 6 5; 2 3 3 6 6 6 6 5 5; 3 3 3 3 6 6 6 4 5; 1 1 5 6 6 6 6 6 6; 1 1 5 5 6 6 6 5 6; 1 6 6 6 5 5 2 1 3; 1 2 2 2 2 4 4 5 3]],
    #                             [[2 3 3 3 3 1 3 3 3; 3 3 3 1 3 1 3 3 3; 6 6 3 2 1 4 2 2 5; 4 4 1 1 4 5 4 4 4; 4 4 4 1 1 1 1 3 5; 4 4 4 1 1 1 3 3 5; 4 6 1 1 4 4 4 4 4; 5 4 4 4 4 3 3 5 5; 5 5 5 5 5 5 6 5 5], [4 4 1 5 5 3 2 4 4; 1 1 1 2 2 3 2 2 2; 1 1 6 2 2 5 3 2 6; 4 2 2 2 2 2 2 2 3; 2 2 2 2 2 2 2 2 2; 2 2 2 2 2 2 2 2 2; 2 4 2 3 3 1 1 2 4; 2 5 5 3 3 3 1 6 2; 1 4 6 4 6 4 4 4 6], [3 2 2 2 1 3 3 2 1; 3 2 3 3 1 2 2 2 1; 5 5 3 3 3 2 2 2 6; 5 5 3 3 1 4 1 1 6; 1 1 2 3 3 5 1 1 4; 1 1 5 4 3 5 1 1 2; 1 1 5 2 4 5 2 2 2; 1 4 5 2 2 2 5 5 6; 1 6 2 2 1 5 5 5 6], [3 3 3 3 3 1 1 1 1; 3 3 5 5 3 4 4 6 6; 1 3 5 6 3 1 4 4 5; 5 5 6 6 6 1 1 4 2; 5 1 6 4 4 5 5 5 4; 5 5 6 4 4 5 5 1 6; 5 5 6 6 6 5 5 6 5; 6 5 6 6 6 5 6 6 5; 4 5 5 6 6 6 2 1 6], [5 4 4 1 3 3 6 6 6; 4 4 4 5 5 6 1 1 6; 4 5 5 5 5 3 1 5 1; 2 5 5 5 3 3 3 6 6; 2 6 5 5 5 3 3 4 6; 2 6 4 4 5 3 4 6 6; 4 4 4 4 4 3 4 6 6; 5 4 4 4 4 4 2 2 2; 5 4 4 4 4 4 2 2 2], [2 6 6 6 6 6 3 1 3; 2 6 6 6 6 5 3 3 3; 3 3 6 6 6 4 3 3 3; 3 3 6 6 4 3 3 3 3; 2 6 6 6 6 6 5 5 3; 1 6 6 6 6 6 5 1 1; 3 6 6 6 1 1 1 1 2; 4 1 1 6 5 1 1 1 6; 2 1 1 1 1 1 1 1 4]],
    #                             [[2 3 3 6 6 6 4 4 4; 6 6 6 6 6 3 5 4 4; 6 6 1 6 1 1 5 2 4; 6 1 1 1 1 1 1 1 4; 1 1 1 1 1 1 2 2 4; 1 1 1 1 1 1 1 2 4; 1 1 1 4 3 1 1 2 2; 1 1 1 2 2 2 1 2 4; 3 3 5 2 6 2 1 4 4], [4 4 6 2 2 2 1 1 1; 4 4 1 4 2 6 6 1 1; 1 2 1 4 2 6 6 5 6; 1 1 2 2 2 6 6 5 6; 4 6 2 2 2 2 6 3 1; 2 6 6 3 3 2 4 3 1; 4 6 6 3 6 5 3 3 1; 4 4 5 5 3 1 3 1 1; 5 6 5 1 1 1 1 1 1], [2 2 2 3 1 3 3 5 5; 2 3 2 5 1 2 4 6 6; 2 2 2 4 1 1 4 6 6; 2 2 3 3 3 3 3 6 6; 2 1 3 3 3 2 1 6 5; 1 1 5 5 5 2 6 6 6; 2 1 5 5 5 1 6 6 6; 2 2 4 4 2 4 5 5 2; 3 5 3 3 3 5 5 5 2], [2 2 1 1 2 6 6 6 1; 2 5 2 3 6 6 6 6 3; 4 1 2 3 6 6 6 6 3; 4 4 4 4 4 4 6 6 6; 6 4 4 4 4 4 4 1 6; 5 2 2 5 3 6 2 1 1; 5 5 5 2 5 2 2 1 1; 5 3 3 5 5 2 2 2 5; 6 6 6 5 5 4 2 2 4], [6 6 6 5 5 5 5 1 1; 1 1 3 4 5 5 5 5 1; 5 4 3 3 5 5 4 4 4; 5 5 5 5 5 4 4 4 4; 5 5 5 5 5 4 4 4 4; 3 5 5 5 6 4 4 4 4; 3 3 3 3 3 3 4 1 2; 3 3 3 3 3 3 4 5 5; 3 3 3 3 3 3 5 5 5], [6 6 4 2 2 3 2 2 6; 3 2 4 2 3 1 5 6 6; 3 3 2 2 3 3 3 3 3; 5 6 2 2 6 3 4 3 3; 3 4 2 5 6 6 6 5 3; 5 5 5 6 6 6 6 3 2; 5 5 5 5 4 2 4 2 4; 5 3 4 4 4 3 4 4 3; 5 4 4 4 4 4 2 3 3]],
    #                             [[4 3 3 3 2 2 2 2 2; 1 3 3 3 5 2 2 2 2; 3 3 3 3 5 5 5 2 2; 3 3 3 3 1 5 2 2 2; 5 6 6 3 1 5 2 2 2; 6 6 1 1 2 2 2 5 1; 6 6 1 1 2 2 2 1 5; 3 6 2 4 6 6 5 1 5; 6 6 3 5 6 3 3 1 5], [3 3 5 5 6 6 1 5 6; 5 1 1 1 5 6 1 6 6; 4 1 1 1 5 5 2 2 1; 1 1 1 1 5 6 6 2 1; 2 1 1 2 2 2 1 4 1; 1 3 5 2 2 2 1 1 1; 5 1 3 4 4 2 3 3 5; 5 5 5 5 3 3 3 3 2; 5 6 5 6 5 5 6 5 6], [4 4 4 4 5 3 3 3 6; 6 6 6 6 6 1 5 2 6; 6 3 2 6 6 1 1 2 2; 6 3 3 5 6 1 3 2 2; 5 3 3 3 3 1 2 2 3; 5 5 5 5 5 1 3 2 2; 5 5 5 4 3 4 4 4 2; 5 5 5 4 4 4 6 4 4; 5 5 4 4 4 4 2 4 4], [3 4 5 2 2 2 2 2 1; 3 5 6 6 2 6 2 2 6; 6 6 6 6 5 6 6 2 6; 5 3 5 3 1 4 4 5 1; 3 1 1 1 4 4 1 1 1; 3 1 1 3 3 4 6 6 4; 3 1 1 3 3 3 6 6 6; 1 1 1 1 1 2 4 6 2; 1 1 1 1 1 6 6 2 2], [2 2 1 6 1 2 2 3 3; 1 4 5 5 5 2 2 2 2; 1 4 5 5 4 2 4 4 4; 3 5 5 5 5 6 4 4 4; 3 6 6 4 5 4 4 4 4; 1 1 1 4 4 4 4 4 4; 1 3 4 4 4 4 4 4 4; 1 3 4 4 4 4 4 4 4; 5 4 4 5 4 4 4 4 4], [1 6 4 5 6 5 5 1 1; 6 4 4 4 3 3 5 5 5; 6 6 6 6 6 6 5 5 1; 6 6 6 6 6 3 5 5 6; 3 5 5 6 6 6 2 2 6; 2 1 2 2 3 6 2 2 3; 2 1 2 2 3 3 3 3 3; 1 1 6 5 3 3 3 3 3; 2 4 1 4 4 3 3 3 3]]]
            

    # E_4_starting_configurations = [ [[1 1 6 6 1 1 6 6 6; 1 3 3 3 1 1 6 6 6; 4 3 3 3 3 4 4 6 6; 4 3 3 3 4 4 4 5 5; 4 3 1 1 1 3 3 4 4; 4 4 1 1 1 1 1 1 4; 1 1 2 2 6 1 1 1 4; 1 1 2 2 6 6 6 1 4; 1 1 2 2 2 1 3 3 1], [2 2 2 2 6 2 2 3 1; 2 2 2 2 2 2 2 2 1; 2 2 2 2 2 2 2 2 5; 1 3 3 2 2 2 2 6 6; 2 3 2 2 2 2 2 6 6; 2 2 2 2 2 2 6 6 6; 2 2 2 2 2 2 4 6 6; 2 2 2 2 2 2 4 3 3; 2 2 2 2 2 2 2 3 3], [6 2 4 3 3 3 1 5 4; 6 3 3 3 2 3 1 1 1; 3 3 3 3 3 3 6 1 1; 3 6 3 3 3 1 1 1 1; 6 6 6 3 3 1 1 1 1; 6 6 5 5 1 1 1 1 1; 3 6 5 2 1 1 1 1 1; 3 4 6 1 1 1 1 1 1; 3 1 1 6 1 1 1 4 2], [4 4 3 3 3 4 4 4 4; 3 4 4 4 4 4 4 4 4; 3 3 1 5 4 4 4 4 4; 3 3 5 3 4 4 4 4 4; 3 3 3 3 4 4 4 4 4; 3 3 3 3 4 4 4 4 6; 4 3 3 4 4 4 4 4 6; 3 3 4 4 4 4 4 4 4; 3 3 3 4 4 4 4 4 3], [5 5 5 4 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 1 5 5 5 5 5 6 5 5; 2 2 5 5 5 5 6 5 5; 2 2 5 5 5 6 6 6 5; 2 2 5 5 5 3 3 5 4; 2 2 5 5 5 3 5 5 4; 2 2 5 5 5 5 5 5 5], [4 6 3 3 3 3 1 6 6; 2 6 3 1 3 4 4 6 6; 3 3 6 6 6 6 6 6 6; 2 2 1 6 6 6 6 5 5; 6 1 5 6 6 6 1 5 5; 6 6 6 6 6 4 1 5 5; 6 5 6 6 4 4 1 1 5; 6 6 6 6 6 1 1 5 5; 6 6 6 6 1 1 1 5 5]],
    #                             [[1 1 1 6 5 1 1 1 1; 1 6 2 1 1 1 1 1 1; 2 2 2 5 5 1 1 1 1; 2 2 2 1 5 1 1 1 1; 3 3 1 1 1 1 1 5 1; 4 5 1 1 1 1 3 3 5; 4 1 1 1 1 3 3 3 3; 5 1 1 1 1 1 3 3 3; 5 6 1 1 1 1 3 3 4], [2 4 6 2 2 5 5 5 4; 2 2 2 2 2 2 2 6 4; 3 6 2 2 2 2 2 4 4; 1 1 6 2 2 2 2 4 4; 1 1 3 2 2 2 2 4 4; 1 1 3 3 1 2 2 4 2; 2 5 5 2 2 2 2 2 5; 2 2 1 2 2 2 2 2 4; 2 2 2 2 2 2 2 2 2], [3 4 6 3 2 3 1 3 3; 2 3 3 3 3 3 3 3 3; 3 3 3 1 3 3 3 3 3; 3 3 3 2 3 3 3 3 6; 3 3 3 3 3 3 1 1 1; 3 3 3 3 3 3 1 1 1; 3 3 6 5 3 6 1 1 6; 3 3 6 5 6 6 1 1 1; 3 3 2 6 6 3 1 1 1], [1 3 6 4 4 4 2 2 3; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 1; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 3; 4 4 4 4 5 4 4 4 4; 2 5 5 6 6 2 5 1 3; 2 5 4 4 3 2 5 5 5], [5 5 5 5 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 6 2 5 5 5 5 5 5 5; 6 2 5 5 5 5 5 5 5; 6 2 2 2 5 5 5 5 5; 6 6 2 5 5 5 5 5 5; 6 6 6 6 4 5 5 5 5; 6 6 6 5 5 5 1 4 4; 6 6 6 6 5 5 1 1 4], [6 6 3 3 3 3 3 6 6; 6 6 6 3 3 3 3 2 6; 6 6 6 6 6 3 3 6 2; 6 6 6 6 6 6 6 6 6; 2 2 6 6 6 6 6 6 6; 2 2 1 6 6 6 6 6 4; 2 6 1 1 6 6 6 4 4; 2 5 5 6 6 6 2 4 4; 6 1 1 2 6 5 5 6 4]],
    #                             [[5 5 5 5 5 5 1 1 5; 5 5 5 5 5 5 5 5 5; 5 6 5 5 5 5 5 5 5; 6 6 6 1 1 5 5 5 5; 6 6 6 3 1 5 5 5 5; 6 6 6 3 1 5 5 5 6; 6 6 6 1 1 5 5 5 1; 1 2 6 1 5 5 5 5 5; 5 5 5 5 5 5 5 5 5], [4 6 4 2 2 2 2 2 2; 4 2 2 2 2 2 2 2 2; 5 5 2 2 2 2 2 6 5; 1 1 2 2 2 2 2 3 5; 1 1 2 2 2 2 2 6 5; 1 4 1 2 2 2 5 1 5; 1 6 1 2 2 2 2 3 6; 1 6 2 2 2 2 2 4 1; 1 2 2 2 2 2 2 4 1], [1 1 1 4 4 4 4 1 1; 2 1 1 6 1 3 4 4 3; 4 2 6 6 3 3 3 3 3; 4 4 3 3 3 3 3 3 3; 4 4 4 4 3 3 3 3 3; 4 4 4 4 3 3 3 2 3; 2 4 4 4 3 3 3 3 3; 2 3 1 6 3 3 3 3 3; 3 3 3 3 3 3 3 3 3], [3 3 3 1 1 1 1 3 3; 1 1 1 1 1 1 1 1 3; 1 1 1 1 1 1 1 1 1; 1 3 1 1 1 1 1 1 1; 1 1 1 1 4 5 4 4 1; 1 1 1 1 4 4 4 4 1; 1 1 1 1 5 3 4 6 6; 1 1 1 6 3 6 4 6 6; 4 4 5 5 6 6 3 6 2], [6 6 2 2 2 2 2 5 2; 3 3 3 2 2 2 2 2 2; 3 3 3 2 1 2 2 2 4; 3 3 3 6 6 6 6 2 2; 3 3 3 5 5 6 6 5 2; 3 3 3 5 5 6 6 5 2; 6 3 3 5 5 5 5 5 2; 6 6 4 5 6 5 5 5 2; 6 6 6 6 4 6 6 5 2], [4 2 2 3 3 3 3 6 6; 4 4 4 4 4 3 3 3 4; 4 4 4 4 4 6 6 6 6; 4 4 4 4 4 5 6 6 6; 4 4 4 4 6 6 6 6 6; 4 4 4 4 6 6 6 6 6; 4 4 4 4 6 4 6 6 6; 4 4 4 4 2 1 2 6 6; 4 4 4 4 6 4 4 4 6]]]


    # T = 1.9952623149688797
    # E_2p5_starting_configurations = [ [[3 3 6 1 2 5 1 1 4; 2 2 2 4 2 2 1 1 4; 2 2 1 6 6 6 1 2 2; 2 1 1 2 4 4 1 1 2; 4 4 4 4 1 6 6 4 1; 4 4 4 4 1 5 6 6 4; 3 3 6 4 3 5 3 3 1; 2 3 3 6 2 2 3 3 2; 4 5 2 2 2 3 1 3 5], [5 5 4 2 4 4 4 4 4; 6 5 6 5 5 4 4 4 3; 6 2 6 2 2 4 2 4 5; 6 1 3 3 1 2 2 2 4; 6 1 1 3 2 2 4 4 2; 5 3 5 1 4 1 1 4 4; 3 3 5 1 4 6 1 4 3; 3 3 5 5 5 6 4 4 4; 2 5 4 2 1 1 2 2 6], [1 5 2 2 2 3 3 3 1; 1 2 2 2 2 3 3 1 6; 4 4 2 5 3 3 4 4 6; 4 4 4 3 3 3 3 3 6; 4 4 2 2 3 3 3 3 6; 4 4 4 4 5 3 5 3 5; 4 4 4 6 6 4 4 6 6; 2 6 6 6 6 4 6 6 3; 2 3 3 1 6 6 1 1 3], [1 2 6 6 3 6 5 4 4; 5 5 5 6 3 5 4 4 6; 5 6 6 6 5 5 6 5 5; 2 6 2 2 4 4 5 5 5; 4 6 2 2 4 5 5 5 5; 2 3 4 6 3 5 5 5 5; 3 5 4 4 5 1 5 5 5; 4 6 5 2 6 6 1 5 5; 2 2 6 6 6 6 6 6 2], [3 1 1 1 1 1 6 6 6; 3 1 1 1 1 4 6 6 6; 1 1 1 1 4 2 3 6 1; 1 1 1 1 1 5 3 5 5; 1 1 1 5 5 2 6 6 5; 1 1 2 5 5 2 3 1 3; 5 1 5 2 2 2 2 1 4; 5 2 2 2 2 2 2 4 6; 5 1 2 3 3 3 4 4 6], [3 5 2 3 5 6 2 2 5; 6 5 6 1 1 2 2 2 4; 5 5 3 1 1 6 5 5 4; 5 3 3 6 6 6 3 5 4; 5 3 1 6 6 6 5 5 3; 5 6 6 1 1 6 5 5 1; 5 1 2 2 3 3 3 3 1; 1 1 1 3 3 3 3 3 4; 1 1 3 3 3 3 3 1 6]],
    #                                 [[4 4 6 5 6 6 6 2 5; 6 4 4 4 5 4 3 3 6; 6 4 4 6 5 3 3 3 4; 5 3 2 3 3 3 4 5 4; 5 3 2 5 1 1 6 6 6; 6 6 5 6 6 6 6 6 6; 1 5 5 5 6 4 4 6 1; 1 1 1 1 1 4 4 4 3; 6 6 1 1 1 1 1 6 5], [6 5 5 3 6 1 3 5 3; 5 5 5 5 6 1 6 5 1; 5 4 5 5 1 1 6 6 6; 5 5 5 1 6 4 3 1 5; 6 6 4 4 2 2 5 1 3; 6 6 2 2 2 4 4 5 3; 5 6 2 2 5 3 4 5 2; 5 2 1 3 3 3 3 5 2; 6 6 3 3 3 3 2 5 2], [4 6 6 6 4 4 4 4 6; 6 6 6 2 4 4 4 6 6; 3 6 6 6 6 2 6 6 1; 5 6 6 6 6 5 1 3 6; 5 5 2 2 3 5 4 2 2; 4 6 2 2 3 1 5 5 2; 3 1 2 1 3 1 1 4 3; 3 1 2 2 1 1 4 4 3; 3 1 2 2 1 4 4 2 2], [3 3 5 5 5 6 5 3 4; 1 3 5 5 5 1 1 6 1; 1 5 5 5 2 3 3 1 2; 4 5 5 5 5 5 5 5 2; 1 2 6 5 4 4 1 5 5; 1 2 1 6 1 4 4 4 5; 6 3 6 6 3 3 3 3 5; 5 3 2 6 3 3 3 3 1; 5 1 2 3 3 3 3 3 1], [2 4 3 1 4 2 4 3 1; 5 4 1 1 4 6 6 6 1; 6 3 1 2 3 3 4 1 1; 3 3 6 4 4 3 6 6 2; 3 3 4 1 5 3 1 6 2; 3 3 4 3 2 2 4 2 1; 2 2 3 3 2 4 2 2 2; 2 1 3 3 2 2 2 2 2; 3 4 4 4 4 5 2 2 4], [1 3 3 2 2 2 1 5 5; 4 2 2 4 4 1 1 1 2; 6 2 2 2 1 1 1 2 5; 6 2 2 2 1 1 1 1 1; 2 2 3 3 6 6 5 1 1; 2 2 3 1 4 5 1 4 4; 4 4 5 4 4 6 1 5 4; 4 5 5 2 4 4 5 2 2; 1 4 5 1 4 4 4 4 2]],
    #                                 [[5 2 1 5 4 4 2 2 2; 4 6 1 5 4 4 3 2 2; 5 5 6 6 4 3 4 6 1; 1 1 1 1 5 3 4 5 3; 1 1 1 1 1 3 5 6 5; 3 2 1 6 1 2 2 2 2; 2 2 2 2 2 3 3 1 5; 3 2 2 2 1 3 3 1 5; 4 4 2 2 4 2 2 1 1], [5 5 1 5 2 2 4 4 3; 6 3 1 2 2 2 3 4 4; 3 3 3 2 2 3 2 2 4; 3 3 3 3 2 3 4 1 4; 3 3 2 2 2 3 3 3 6; 3 1 2 5 2 6 5 5 6; 5 5 5 5 5 6 6 6 6; 5 3 2 5 2 6 5 5 4; 3 3 2 3 1 1 6 6 3], [2 3 6 6 2 3 5 5 3; 2 5 4 4 2 6 6 5 5; 2 2 3 1 1 1 5 5 5; 1 1 4 1 1 1 1 6 5; 4 4 3 3 3 5 5 6 6; 4 4 3 2 3 5 5 4 1; 6 4 6 6 6 5 5 4 5; 4 4 1 6 6 4 5 5 5; 4 3 4 4 3 1 3 5 4], [1 1 3 6 5 5 6 6 6; 1 1 1 1 4 4 5 6 4; 6 2 5 4 3 3 1 6 4; 6 2 2 4 5 3 5 5 1; 2 5 5 4 4 4 1 1 5; 2 5 5 4 4 4 4 3 3; 2 6 2 5 4 5 4 5 3; 2 2 3 3 2 6 2 3 3; 5 5 4 4 1 2 2 1 2], [1 3 5 4 4 4 4 1 4; 1 3 4 4 4 4 4 4 3; 1 5 4 4 6 6 1 3 3; 5 6 6 6 6 4 1 2 4; 6 5 6 6 5 1 1 5 5; 6 6 6 5 2 2 1 5 2; 5 3 6 6 6 3 4 4 1; 1 6 6 6 6 2 4 2 2; 1 1 3 6 6 1 4 6 6], [6 6 3 3 3 6 3 6 6; 6 6 3 3 5 1 4 4 2; 6 6 3 3 4 4 1 6 1; 5 3 4 6 6 2 2 3 2; 2 3 4 6 6 4 2 1 1; 5 3 2 5 5 1 1 1 1; 4 2 2 2 3 6 1 1 1; 4 1 1 5 3 1 1 1 2; 2 3 1 5 3 6 6 6 5]],
    #                                 [[5 5 1 1 1 1 2 3 4; 5 3 5 6 6 3 1 5 1; 6 5 5 5 6 4 4 1 2; 2 2 2 1 2 2 4 2 6; 2 2 2 1 1 1 1 1 1; 1 5 6 5 1 1 1 2 3; 1 1 1 1 1 5 3 4 3; 1 1 1 5 5 3 3 4 4; 2 1 1 5 5 4 4 4 4], [4 2 3 3 5 5 5 5 5; 4 5 3 2 6 6 5 5 5; 5 5 5 5 2 2 5 4 1; 5 5 5 2 2 4 1 1 1; 4 4 2 2 2 4 5 6 2; 4 4 2 2 2 4 4 4 6; 2 4 2 2 2 2 2 3 6; 2 2 2 2 3 3 5 5 6; 2 2 5 2 5 5 1 5 3], [2 2 5 5 6 3 4 1 1; 2 2 2 5 6 6 6 6 6; 2 2 3 3 3 6 6 6 6; 3 3 3 3 3 6 6 6 2; 4 4 4 3 3 5 6 2 3; 4 1 1 6 5 5 3 1 1; 6 3 3 4 4 4 4 6 3; 6 6 6 6 2 4 6 6 5; 6 6 6 6 4 1 4 4 1], [6 4 4 4 3 3 3 3 3; 4 4 2 3 3 3 3 3 3; 4 4 6 3 3 3 5 5 5; 2 4 2 3 6 3 5 5 5; 3 5 5 5 4 4 5 5 5; 3 2 2 4 4 5 5 5 5; 5 4 2 6 5 4 4 3 5; 6 1 2 1 2 4 2 4 6; 6 1 2 4 2 2 6 5 5], [1 1 3 3 1 1 1 3 3; 1 1 1 3 3 1 1 3 3; 1 4 1 1 3 6 6 5 4; 6 6 1 1 3 3 3 5 6; 6 3 6 6 5 3 3 4 1; 6 6 6 6 1 4 4 4 4; 3 6 1 1 1 4 4 4 4; 3 3 6 4 1 1 3 2 2; 3 3 4 4 6 3 3 3 1], [5 1 2 2 3 6 6 4 4; 5 1 1 2 1 1 1 2 2; 3 5 3 3 1 6 6 2 5; 4 5 5 5 6 6 6 4 6; 6 1 4 6 6 5 6 4 4; 5 6 2 2 4 1 5 2 1; 6 6 1 3 4 1 2 2 2; 6 6 3 3 5 1 4 4 2; 6 6 1 2 2 2 2 4 2]],
    #                                 [[3 3 4 4 3 2 2 6 1; 3 3 4 1 4 6 4 5 3; 6 3 4 1 5 4 5 6 5; 4 3 6 1 1 3 3 3 3; 4 3 4 1 1 1 3 6 3; 6 3 3 1 2 2 2 4 3; 5 5 4 4 1 6 2 2 3; 4 1 4 4 4 6 1 1 3; 5 2 4 4 4 5 1 3 3], [6 6 1 1 1 4 3 3 6; 6 6 6 4 1 4 1 3 4; 6 6 6 4 5 6 3 6 4; 6 3 4 4 1 5 3 6 4; 6 5 4 3 2 2 3 6 4; 1 1 2 3 6 6 5 1 5; 2 2 2 3 3 6 5 5 5; 2 6 6 6 6 6 1 5 5; 6 6 6 6 6 6 6 1 5], [5 5 6 2 2 1 1 4 2; 5 5 1 2 2 2 2 4 4; 2 5 6 2 2 1 1 4 4; 4 5 5 2 2 1 1 4 2; 5 6 1 6 3 3 6 1 6; 3 5 2 6 3 6 1 1 6; 5 5 3 6 6 6 6 5 5; 3 3 3 4 4 6 3 5 1; 3 4 5 5 5 6 4 5 4], [4 6 3 5 1 1 4 5 5; 2 3 4 2 3 5 4 2 5; 6 2 2 2 3 3 1 3 1; 5 2 2 3 4 3 1 1 2; 3 2 2 3 4 4 1 3 3; 5 5 5 1 4 4 6 6 2; 2 6 1 1 4 4 6 6 1; 4 2 1 1 5 2 2 1 1; 4 4 4 1 1 6 6 1 6], [2 2 2 2 2 2 3 1 4; 2 2 2 2 2 2 6 1 6; 2 2 5 5 2 2 1 1 1; 2 2 5 5 5 2 3 3 1; 2 5 5 5 5 5 6 5 5; 3 3 5 5 5 5 3 3 3; 3 3 2 5 5 5 3 3 3; 1 6 1 3 4 5 3 4 3; 1 6 1 3 5 4 3 1 3], [2 2 2 3 6 1 1 1 1; 4 2 2 4 3 5 1 6 6; 4 4 4 6 4 2 5 5 3; 6 1 1 4 2 2 1 1 1; 4 1 6 6 6 4 1 1 1; 4 4 4 6 6 4 4 5 5; 6 4 3 3 2 4 4 5 5; 2 4 3 6 2 5 5 4 5; 2 2 2 3 2 5 5 5 1]]]


    # T = T=1.9952623149688797
    starting_configurations = [[[4 6 6 5 1 4 4 4 1; 6 6 1 1 1 5 4 4 5; 6 3 6 6 5 5 5 5 5; 6 6 6 6 5 5 5 5 5; 6 5 5 6 1 2 4 4 5; 6 2 2 2 4 1 3 4 4; 5 5 1 1 3 3 3 3 4; 5 5 1 2 2 3 3 2 6; 5 6 1 3 1 1 6 1 1], [2 5 5 6 2 2 2 2 2; 2 3 2 6 2 2 2 6 1; 6 6 2 3 2 2 2 4 1; 6 6 6 6 2 2 1 4 1; 1 1 6 6 2 1 6 3 3; 3 1 6 3 2 1 6 4 3; 6 2 2 3 1 1 6 6 6; 2 2 2 3 1 1 1 6 5; 5 1 4 3 4 5 1 6 6], [2 2 2 2 2 1 5 3 3; 2 2 2 2 3 1 5 5 2; 2 5 5 6 6 2 4 4 2; 6 5 5 4 3 6 4 4 4; 6 6 3 3 3 3 4 4 4; 6 6 1 3 3 1 1 1 2; 6 6 6 3 3 3 3 1 3; 1 1 3 3 3 3 3 3 3; 1 1 1 1 3 3 3 3 3], [3 4 4 2 2 2 2 5 5; 5 2 4 2 2 2 2 3 5; 5 4 4 2 2 2 5 3 3; 5 4 4 3 1 1 5 5 5; 3 3 5 6 4 1 3 5 5; 3 3 5 6 4 5 3 5 5; 3 5 5 5 4 4 2 5 5; 3 5 3 3 6 4 6 4 1; 5 5 3 1 6 4 4 4 4], [6 6 1 1 5 5 3 3 3; 4 4 1 5 5 5 2 1 2; 1 1 1 5 5 5 3 6 2; 4 4 1 5 5 5 1 1 2; 3 4 1 5 5 1 1 1 1; 3 4 4 4 2 2 2 1 3; 1 1 4 4 2 2 1 1 3; 1 1 4 2 2 3 3 3 3; 1 4 4 4 5 4 3 3 2], [4 4 4 4 4 5 5 6 6; 2 5 4 6 6 1 6 1 6; 2 5 3 3 4 4 1 5 4; 2 3 1 3 4 4 4 6 6; 4 4 1 6 6 4 6 6 6; 1 5 6 2 5 4 6 6 6; 5 6 6 4 2 2 4 6 1; 3 6 4 6 5 2 2 4 1; 4 4 6 1 2 2 2 4 6]],
                                [[5 5 5 2 2 6 6 4 4; 5 5 3 5 5 5 1 3 3; 5 5 3 5 4 4 6 3 2; 5 5 1 5 5 4 2 3 3; 5 5 2 2 1 1 2 2 2; 5 5 2 2 5 2 2 2 2; 4 6 2 2 3 1 1 6 1; 4 2 2 2 2 3 5 6 1; 1 2 2 4 2 3 3 3 3], [3 4 2 2 6 1 4 2 2; 6 2 2 2 2 1 5 2 1; 6 5 5 1 1 1 1 1 1; 1 1 5 5 1 1 1 1 1; 3 4 5 4 2 5 4 4 6; 3 5 5 1 6 6 4 4 4; 3 3 5 5 6 6 4 4 3; 3 3 3 5 5 1 4 5 3; 3 3 3 3 5 1 3 6 5], [2 2 1 5 5 3 1 1 1; 5 2 2 2 3 3 1 1 3; 6 5 4 4 1 1 1 1 3; 6 4 4 6 3 2 3 4 3; 4 1 3 3 3 1 1 2 3; 4 3 3 3 3 1 1 6 1; 6 6 4 3 3 4 1 1 1; 1 3 1 1 1 1 5 1 1; 6 3 3 1 1 5 6 6 3], [2 1 1 1 1 2 2 5 5; 3 5 1 5 1 1 2 4 5; 4 6 5 5 5 2 2 2 5; 6 6 6 5 4 1 1 2 1; 4 4 4 4 4 1 1 3 3; 5 6 6 4 4 3 3 3 3; 5 4 4 4 4 3 3 3 1; 5 6 3 4 4 3 3 6 6; 5 5 5 4 4 4 4 6 6], [4 4 4 6 5 6 5 6 4; 4 4 6 1 1 2 4 4 4; 2 2 2 5 2 2 3 4 4; 5 3 3 5 5 4 4 4 2; 1 5 2 2 5 2 5 3 1; 2 2 2 3 3 3 5 5 2; 2 2 2 3 3 3 3 3 3; 1 1 4 3 3 4 4 6 1; 1 4 4 3 3 4 1 6 6], [1 2 6 6 4 5 6 6 6; 2 3 6 6 6 6 6 4 5; 5 6 6 6 6 6 6 5 5; 5 6 6 6 6 4 4 4 4; 6 6 6 6 6 6 6 6 2; 6 6 6 6 2 2 6 4 4; 6 1 5 5 5 2 6 4 4; 2 1 2 2 6 6 5 5 4; 2 2 2 2 6 6 2 2 4]],
                                [[3 3 3 3 1 3 3 3 3; 3 3 4 3 3 3 3 3 3; 3 2 4 1 1 1 1 3 3; 6 6 6 6 1 1 1 2 4; 2 6 6 6 1 1 5 4 4; 6 6 3 4 1 1 1 1 4; 3 2 5 2 2 2 4 4 4; 3 2 5 2 2 2 5 4 4; 6 5 6 2 6 2 4 4 1], [5 1 4 4 4 4 2 2 2; 5 3 6 1 3 4 2 2 2; 1 5 6 3 3 2 2 2 2; 3 3 6 3 3 2 2 3 3; 1 3 3 3 2 2 2 2 2; 1 2 6 3 2 2 4 4 2; 1 6 3 3 3 6 1 1 1; 6 6 3 3 2 3 6 6 6; 6 6 4 2 2 5 2 5 6], [5 5 5 4 4 1 1 1 1; 2 5 5 5 1 1 1 4 1; 2 4 5 5 2 1 1 1 1; 4 4 3 1 4 1 1 1 1; 4 4 4 5 3 1 2 2 1; 6 1 2 2 2 2 2 2 1; 6 3 5 4 1 1 6 6 2; 6 4 4 4 6 6 6 6 5; 3 4 6 6 6 1 1 3 5], [1 1 1 1 3 5 5 1 4; 3 1 1 1 1 1 1 1 4; 3 3 2 4 1 1 1 2 2; 3 6 3 4 6 6 6 6 2; 5 6 1 4 4 6 6 1 3; 5 2 2 6 4 6 6 6 1; 6 1 2 6 5 6 6 6 4; 1 2 2 2 5 4 3 5 4; 4 2 3 3 1 2 5 5 5], [4 4 4 2 6 6 4 3 4; 1 5 4 4 4 6 4 4 4; 1 1 3 3 3 2 4 4 4; 4 1 5 4 3 3 4 4 3; 5 5 5 5 5 2 4 6 3; 3 5 5 5 3 3 3 3 1; 3 3 4 4 4 3 3 3 6; 6 1 6 2 4 6 1 3 5; 1 2 2 2 2 6 6 2 2], [2 2 2 6 6 5 5 1 3; 4 2 2 5 3 5 5 1 2; 5 2 2 5 4 4 3 5 5; 5 5 5 5 4 4 4 5 5; 3 1 6 5 6 6 6 5 5; 5 4 4 5 5 5 5 5 5; 5 4 6 5 5 5 5 5 5; 5 6 6 5 5 3 5 5 6; 6 6 6 6 5 4 6 6 2]]]

    for (index,starting_configuration) in pairs(starting_configurations)
        cube = RubiksCube(9)
        cube.configuration = starting_configuration

        time_experiment(cube, "funnel_test_E=$(E)_T=$(T)_trial_$index", 9, [p_swap], T, N_t; normalization=normalization)

    end

end