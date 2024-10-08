using DelimitedFiles
using SharedArrays
using Distributed


include("../probes/relaxed_anneal.jl")
include("../tools/neighbour_graphs_plotter.jl")



@inbounds @fastmath function connections_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, relaxation_iterations::Int64=Int(0), connections_to_measure::Union{String,Nothing}="slice", connections_per_configuration_sample_size::Int64=0, neighbour_order_to_measure_to::Int64=1, cumulative_connections_measurement::Bool=false, average_sample_size_per_temperature::Int64=100, inherent_disorder::Bool=false, T_1::Float64=10.0, T_0::Float64=0.09, T_swap::Float64=3.0, initial_cube_configuration=nothing, parallel_anneals::Int64=1, collect_energy_saddle_index_densities::Bool=true, bin_energy_connections::Bool=false, disorder_average::Bool=false)

    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; sample_temperatures]

    sort!(temperature_vector, rev=true)
    sample_temperatures = sort!(sample_temperatures)


    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    if isnothing(initial_cube_configuration)
        if inherent_disorder
            facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
            shuffle!(facelets)
            new_faces = reshape(facelets, 6, L, L)
            for i in 1:6
                cube.configuration[i][:,:] .= new_faces[i,:,:]
            end
        end
    else
        println("Using initial cube configuration")
        cube.configuration = initial_cube_configuration
    end

    # Evaluate one or many relaxed_anneals depending on if using parallel processing or not
    if parallel_anneals > 1

        samples_per_temperature_per_anneal = Int(ceil(average_sample_size_per_temperature/parallel_anneals))

        Z = configuration_network_degree(cube.L, connections_to_measure=="swap") # Note this TESTS a bool
        number_of_connections = connections_per_configuration_sample_size==0 ? Z*(Z-1)^(neighbour_order_to_measure_to-1) : connections_per_configuration_sample_size

        # Calculate the sizes for each SharedArray
        size_energy_connections = length(sample_temperatures) * number_of_connections * samples_per_temperature_per_anneal * parallel_anneals
        size_energy_saddle_index_density = length(sample_temperatures) * samples_per_temperature_per_anneal * parallel_anneals
        size_energy_minima = length(sample_temperatures) * samples_per_temperature_per_anneal * parallel_anneals

        println("Size of energy connections: $size_energy_connections")

        # Initialize SharedArrays with the correct sizes
        energy_connections_data = cumulative_connections_measurement ? SharedArray{Float64,2}((size_energy_connections, neighbour_order_to_measure_to+1))        : SharedArray{Tuple{Float64,Float64},1}((size_energy_connections,))
        energy_saddle_index_density_tuples = SharedArray{Tuple{Float64,Float64},1}((size_energy_saddle_index_density,))
        energy_minima_tuples = SharedArray{Tuple{Float64,Bool},1}((size_energy_minima,))

        @sync @distributed for trial in 1:parallel_anneals

            # Initialize the cube with different inherent disorder per trial if disorder_average
            if disorder_average
                facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
                shuffle!(facelets)
                new_faces = reshape(facelets, 6, L, L)
                for i in 1:6
                    cube.configuration[i][:,:] .= new_faces[i,:,:]
                end
            end

            printstyled("Trial: $trial \n", color=:light_blue)

            if relaxation_iterations == 0
                temperature_vector, _, _, _, _, _, _, _, _, trial_energy_connections_data, trial_energy_saddle_index_densities_tuples, trial_energy_minima_tuples = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, sample_temperatures=sample_temperatures, connections_to_measure=connections_to_measure, connections_per_configuration_sample_size=connections_per_configuration_sample_size, bin_energy_connections=bin_energy_connections, cumulative_connections_measurement=cumulative_connections_measurement, collect_energy_saddle_index_densities=collect_energy_saddle_index_densities, average_sample_size_per_temperature=samples_per_temperature_per_anneal, neighbour_order_to_measure_to=neighbour_order_to_measure_to)
            else
                relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
                temperature_vector, _, _, _, _, _, _, _, _, trial_energy_connections_data, trial_energy_saddle_index_densities_tuples, trial_energy_minima_tuples = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector,  sample_temperatures=sample_temperatures, connections_to_measure=connections_to_measure, connections_per_configuration_sample_size=connections_per_configuration_sample_size, bin_energy_connections=bin_energy_connections, cumulative_connections_measurement=cumulative_connections_measurement, collect_energy_saddle_index_densities=collect_energy_saddle_index_densities, average_sample_size_per_temperature=samples_per_temperature_per_anneal, neighbour_order_to_measure_to=neighbour_order_to_measure_to)
            end

            if !bin_energy_connections
                if !cumulative_connections_measurement
                    energy_connections_data[(trial-1)*length(trial_energy_connections_data)+1:trial*length(trial_energy_connections_data)] .= trial_energy_connections_data[:]
                else
                    energy_connections_data[(trial-1)*size(trial_energy_connections_data,1)+1:trial*size(trial_energy_connections_data,1),:] .= trial_energy_connections_data
                end
            end

            if collect_energy_saddle_index_densities
                energy_saddle_index_density_tuples[(trial-1)*length(trial_energy_saddle_index_densities_tuples)+1:trial*length(trial_energy_saddle_index_densities_tuples)] .= trial_energy_saddle_index_densities_tuples[:]
                energy_minima_tuples[(trial-1)*length(trial_energy_minima_tuples)+1:trial*length(trial_energy_minima_tuples)] .= trial_energy_minima_tuples[:]
            end
        end
    else

        if relaxation_iterations == 0
            temperature_vector, _, _, _, _, _, _, _, _, trial_energy_connections_data, trial_energy_saddle_index_densities_tuples, trial_energy_minima_tuples = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, sample_temperatures=sample_temperatures, connections_to_measure=connections_to_measure, connections_per_configuration_sample_size=connections_per_configuration_sample_size, bin_energy_connections=bin_energy_connections, cumulative_connections_measurement=cumulative_connections_measurement, collect_energy_saddle_index_densities=collect_energy_saddle_index_densities, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_order_to_measure_to=neighbour_order_to_measure_to)
        else
            relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
            temperature_vector, _, _, _, _, _, _, _, _, trial_energy_connections_data, trial_energy_saddle_index_densities_tuples, trial_energy_minima_tuples = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, sample_temperatures=sample_temperatures, connections_to_measure=connections_to_measure, connections_per_configuration_sample_size=connections_per_configuration_sample_size, bin_energy_connections=bin_energy_connections, cumulative_connections_measurement=cumulative_connections_measurement, collect_energy_saddle_index_densities=collect_energy_saddle_index_densities, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_order_to_measure_to=neighbour_order_to_measure_to)
        end

        energy_connections_data = trial_energy_connections_data
        energy_saddle_index_density_tuples = trial_energy_saddle_index_densities_tuples
        energy_minima_tuples = trial_energy_minima_tuples

    end

    ##### ----- ANALYSE RESULTS ----- 

    printstyled("Analyzing results \n", color=:light_blue)

    if !cumulative_connections_measurement
        # Convert the array of tuples to a multi-column matrix
        energy_connections_data_matrix = hcat([x[1] for x in energy_connections_data], [x[2] for x in energy_connections_data])
    else
        energy_connections_data_matrix = energy_connections_data
    end
    # Make energies with respect to solved energy
    energy_connections_data_matrix[:] .= energy_connections_data_matrix[:] .- solved_configuration_energy(cube)


    # Convert the array of tuples to a 2-column matrix
    energy_saddle_index_densities_data_matrix = hcat([x[1] for x in energy_saddle_index_density_tuples], [x[2] for x in energy_saddle_index_density_tuples])
    # Make energies with respect to solved energy
    energy_saddle_index_densities_data_matrix[:,1] .= energy_saddle_index_densities_data_matrix[:,1] .- solved_configuration_energy(cube)


    # Convert the array of tuples to a 2-column matrix
    energy_minima_data_matrix = hcat([x[1] for x in energy_minima_tuples], [x[2] for x in energy_minima_tuples])
    # Make energies with respect to solved energy
    energy_minima_data_matrix[:,1] .= energy_minima_data_matrix[:,1] .- solved_configuration_energy(cube)



    ### --- Save the results ---
    try
        data_simulation_name_to_use = "$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"

        # Energy Connection Results
        touch(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use))

        open(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Average Sample Size = $(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(connections_per_configuration_sample_size), Collecting Swap Move Neighbours = $(connections_to_measure) \n")
            
            if !cumulative_connections_measurement
                write(simulation_file, "Initial Energy E$(neighbour_order_to_measure_to-1), Neighbour Energy E$(neighbour_order_to_measure_to) \n")
                
                for index in 1:size(energy_connections_data_matrix,1)
                    write(simulation_file, "$(energy_connections_data_matrix[index,1]),$(energy_connections_data_matrix[index,2]) \n")
                end
            else
                # Write title line
                title = join(["E$i" for i in 0:neighbour_order_to_measure_to], ",")
                write(simulation_file, "$title\n")
            
                # Write data lines
                for index in 1:size(energy_connections_data_matrix, 1)
                    data_line = join(energy_connections_data_matrix[index, :], ",")
                    write(simulation_file, "$data_line\n")
                end
            end
        end

        data_simulation_name_to_use = "$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_saddle_index_densities"

        # Energy Saddle Index Density Results
        touch(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use))

        open(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Average Sample Size = $(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(connections_per_configuration_sample_size), Collecting Swap Move Neighbours = $(connections_to_measure) \n")
            write(simulation_file, "So saddle index density k = Saddle Index / Neighbours Per Configuration Sample Size \n")
            write(simulation_file, "Energy E$(neighbour_order_to_measure_to-1), Saddle Index Density \n")
            
            for index in 1:size(energy_saddle_index_densities_data_matrix,1)
                write(simulation_file, "$(energy_saddle_index_densities_data_matrix[index,1]),$(energy_saddle_index_densities_data_matrix[index,2]) \n")
            end
        end

        data_simulation_name_to_use = "$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_minima"
        
        # Energy Minima Results
        touch(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use))

        open(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Average Sample Size = $(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(connections_per_configuration_sample_size), Collecting Swap Move Neighbours = $(connections_to_measure) \n")
            write(simulation_file, "Energy E$(neighbour_order_to_measure_to-1), Is Minimum \n")
            
            for index in 1:size(energy_minima_data_matrix,1)
                write(simulation_file, "$(energy_minima_data_matrix[index,1]),$(energy_minima_data_matrix[index,2]) \n")
            end
        end


    catch ex

        println("Cannot save results to file")
        showerror(stdout, ex)

    end

end
