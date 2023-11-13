using Plots
using DelimitedFiles
using LaTeXStrings

using SharedArrays


include("../probes/relaxed_anneal.jl")
include("../tools/neighbour_initial_and_final_energies_graphs_plotter.jl")


# Parallel Case
# using Distributed
# addprocs(Sys.CPU_THREADS-1)
# using SharedArrays
# @everywhere include("../probes/relaxed_anneal.jl")
# @everywhere include("../tools/neighbour_initial_and_final_energies_graphs_plotter.jl")


@inbounds @fastmath function neighbour_initial_and_final_energies_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, collecting_swap_move_neighbours::Bool=false, neighbours_per_configuration_sample_size::Int64=0, average_sample_size_per_temperature::Int64=100, inherent_disorder::Bool=false, neighbour_moves_away::Int64=1)

    sample_temperatures = sort!(sample_temperatures)

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 3.0
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; sample_temperatures]

    sort!(temperature_vector, rev=!temperature_increasing_anneal)
    sample_temperatures = sort!(sample_temperatures)


    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    if inherent_disorder
        facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
        shuffle!(facelets)
        new_faces = reshape(facelets, 6, L, L)
        for i in 1:6
            cube.configuration[i][:,:] .= new_faces[i,:,:]
        end
    end

            if relaxation_iterations == 0
                temperature_vector, _, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=samples_per_temperature_per_anneal, neighbour_moves_away=neighbour_moves_away)
            else
                relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
                temperature_vector, _, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=samples_per_temperature_per_anneal, neighbour_moves_away=neighbour_moves_away)
            end

            neighbour_initial_and_final_energies[(trial-1)*length(neighbour_initial_and_final_energies_by_temperature)+1:trial*length(neighbour_initial_and_final_energies_by_temperature)] .= neighbour_initial_and_final_energies_by_temperature[:]
        end

    else

        if relaxation_iterations == 0
            temperature_vector, _, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_moves_away=neighbour_moves_away, temperature_increasing_anneal=temperature_increasing_anneal)
        else
            relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
            temperature_vector, _, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_moves_away=neighbour_moves_away, temperature_increasing_anneal=temperature_increasing_anneal)
        end

        # First combine (E_1, E_2) samples over all temperatures
        neighbour_initial_and_final_energies = neighbour_initial_and_final_energies_by_temperature[:]

    end

    ##### ----- ANALYSE RESULTS -----

    printstyled("Analyzing results \n", color=:light_blue)

    # First combine (E_1, E_2) samples over all temperatures
    neighbour_initial_and_final_energies = neighbour_initial_and_final_energies_by_temperature[:]

    # Convert the array of tuples to a 2-column matrix
    data_matrix = hcat([x[1] for x in neighbour_initial_and_final_energies], [x[2] for x in neighbour_initial_and_final_energies])

    # Make energies with respect to solved energy
    data_matrix[:] .= data_matrix[:] .- solved_configuration_energy(cube)

    data_simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)"

    ### --- Save the results ---
    try
        touch(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use*"_raw_energy_connections"))

        open(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Average Sample Size = $(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(neighbours_per_configuration_sample_size), Collecting Swap Move Neighbours = $(collecting_swap_move_neighbours) \n")
            write(simulation_file, "Initial Energy E1, Neighbour Energy E2 \n")
            
            for index in 1:size(data_matrix,1)
                write(simulation_file, "$(data_matrix[index,1]),$(data_matrix[index,2]) \n")
            end
        end

    catch ex

        println("Cannot save results to file")
        showerror(stdout, ex)

    end



    ### --- Create and Display the Graphs ---
    connectivity = collecting_swap_move_neighbours ? "Swap" : "Slice"
    neighbour_initial_and_final_energies_graph_plotter(data_simulation_name_to_use, neighbour_moves_away, connectivity)

end
