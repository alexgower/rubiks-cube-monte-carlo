
include("../core/rubiks_cube.jl")
include("../core/monte_carlo.jl")





@inbounds @fastmath function relaxed_anneal!(cube::RubiksCube, temperature_vector::Vector{Float64}; swap_move_probability::Float64=0.0, T_swap::Float64=0.0, relaxation_iterations_vector=nothing, average_sample_size_per_temperature::Int64=100, verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false, mixing_p_swap::Float64=0.0, neighbour_sample_temperatures::Vector{Float64}=empty([0.0]), collecting_swap_move_neighbours::Bool=false, neighbours_per_configuration_sample_size::Int64=0, energy_histogram_sample_temperatures::Vector{Float64}=empty([0.0]), collect_minimum_neighbour_energy_info_only::Bool=false, collect_initial_and_final_energies::Bool=false, neighbour_moves_away::Int64=1)

    # Notes ---

    # Runs temperature anneal on Rubik's Cube using the Metropolis+Swap algorithm at each temperature to
    # relax and probe properties of the cube.



    # Validation and Initial Set-Up --- 
    collecting_neighbour_energy_info = !isempty(neighbour_sample_temperatures)
    creating_energy_histogram = !isempty(energy_histogram_sample_temperatures)


    # Make sure relaxation_iterations_vector (if provided) has same number of elements as temperature_vector
    if !isnothing(relaxation_iterations_vector) && length(relaxation_iterations_vector) != length(temperature_vector)
        throw(ArgumentError("relaxation_iterations_vector must be same length as temperature_vector"))
    end

    # Create default relaxation_iterations_vector if not provided
    if isnothing(relaxation_iterations_vector)
        T1 = temperature_vector[1]
        T0 = temperature_vector[end]

        tau_scaling_factor = 15

        # Create largest relaxation iterations value (for low temperatures) that scales with cube size
        tau_0 = tau_scaling_factor * 6 * cube.L^2

        # Create smallest relaxation iterations (for high temperatures) value that scales with order of generators of cube group
        # Number of generators of size-L Rubik's cube is number of faces * number of layers per face * number of
        # rotation orientations per layer = 6 * ceil((L-1)/2) * 2
        tau_1 = 6 * ceil((cube.L-1)/2) * 2

        relaxation_iterations_vector = [tau_1 * (tau_0/tau_1)^((T1-T)/(T1-T0)) for T in temperature_vector]
    end

    # Variables for later
    tau_0 = relaxation_iterations_vector[1]
    tau_1 = relaxation_iterations_vector[end]


    # Mixing Stage ---

    # Randomise Rubik's cube initially by running Metropolis+Swap algorithm at infinite temperature
    # i.e. beta = 0, swap_move_probability = mixing_p_swap, and by default have maximum_iterations = 10*tau_1
    # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with t=0
    # configuration) has dropped to e^(-10) (i.e. 10 relaxation times) or 10*tau_1 (which should be a reasonable upper bound to this) 
    # iterations have been reached
    run_metropolis_swap_algorithm!(cube, 0.0; swap_move_probability=mixing_p_swap, maximum_iterations=Int(ceil(10*tau_1)), verbose=false, configuration_correlation_convergence_criteria=exp(-10))

    if verbose_annealing
        printstyled("New Cube With P_swap = $swap_move_probability below T_swap = $T_swap \n"; color=:blue)
        println("Mixed cube")
        println("Cube Energy/Infinite Temperature Energy: $(energy(cube)/infinite_temperature_energy(cube))")
    end



    # Annealing Stage ----

    # Create arrays to store parameters for each temeprature
    E_average_by_temperature = zeros(length(temperature_vector))
    E_squared_average_by_temperature = zeros(length(temperature_vector))
    measured_relaxation_iterations_by_temperature = zeros(length(temperature_vector))
    accepted_candidates_by_temperature = zeros(length(temperature_vector))
    final_configuration_correlation_function_by_temperature = zeros(length(temperature_vector))

    # If collecting neighbouring energy deltas then create array to store them
    if collecting_neighbour_energy_info

        # If neighbour_per_configuration_sample_size=0 (default) then just collect all neighbours, otherwise collect a random sample of neighbours
        if neighbours_per_configuration_sample_size==0
            # For Z = degree of network, this is equal to Z*(Z-1)^(neighbour_moves_away-1) since we have Z choices for the first move, and then Z-1 choices for each subsequent move
            Z = configuration_network_degree(cube.L, collecting_swap_move_neighbours)
            number_of_neighbours = Z*(Z-1)^(neighbour_moves_away-1)
        else
            number_of_neighbours = neighbours_per_configuration_sample_size
        end

        neighbour_info_values_per_temperature = collect_minimum_neighbour_energy_info_only ? 1 : number_of_neighbours*average_sample_size_per_temperature

        if collect_initial_and_final_energies
            neighbour_information_by_temperature = Array{Tuple{Float64, Float64}}(undef, length(neighbour_sample_temperatures), neighbour_info_values_per_temperature)
        else
            neighbour_information_by_temperature = zeros(length(neighbour_sample_temperatures),neighbour_info_values_per_temperature)
        end

    end

    # If creating energy histogram then create array to store it (with length just equal to average_sample_size)
    if creating_energy_histogram
        energy_samples_by_temperature = zeros(length(energy_histogram_sample_temperatures),average_sample_size_per_temperature)
    end

    
    # Cool Rubik's cube from T_1 to T_0 by temperatures described in the temperature vector
    for (temperature_index, T) in pairs(temperature_vector)
        try
        beta = 1/T



        # Relaxation Stage ---

        # Run Metropolis+Swap algorithm to relax Rubik's cube at this temperature
        # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with t=0
        # configuration) has dropped to e^(-2) (i.e. 2 relaxation times) or 2*relation_iterations have been reached

        # Only allow swap moves below a certain temperature T_swap 
        if T <= T_swap && swap_move_probability!=0.0
            swap_move_probability_at_this_temperature = swap_move_probability
        else
            swap_move_probability_at_this_temperature = 0.0
        end

        relaxation_converged, final_configuration_correlation_function, final_iteration_number, final_accepted_candidates_number = run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability_at_this_temperature, maximum_iterations=Int(ceil(2*relaxation_iterations_vector[temperature_index])), verbose=verbose_metropolis_swap, configuration_correlation_convergence_criteria=exp(-2))



        # Measurement Stage ---

         # If creating energy histogram at this temperature then measure and store them
        if creating_energy_histogram 
            if insorted(T, energy_histogram_sample_temperatures) && verbose_annealing
                printstyled("!!! Creating energy histogram at T = $T\n"; color=:red)
            end

            # If completed final energy histogram then break out of for loop
            if all(energy_histogram_sample_temperatures .> T)
                break
            end
        end

        # If collecting neighbouring energy deltas at this temperature then measure and store them
        if collecting_neighbour_energy_info
            if insorted(T, neighbour_sample_temperatures) && verbose_annealing
                printstyled("!!! Collecting neighbour energy deltas at T = $T\n"; color=:red)

                if collect_minimum_neighbour_energy_info_only
                    printstyled("!!! Collecting minimum neighbour energy delta only"; color=:red)
                end
            end
            
            # If completed final neighbour energy deltas then break out of for loop
            if all(neighbour_sample_temperatures .> T)
                break
            end
        end

        
        # Do measurements ---

        # Calculate <E> and <E^2> at this temperature but only take measurements every relation_iterations steps to ensure statistical independence
        E_running_total = 0.0
        E_squared_running_total = 0.0


        for sample_index in 1:average_sample_size_per_temperature

            # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with
            # t=0 configuration) has dropped to e^(-1) (i.e. 1 relaxation time) or tau(T) (which should be a reasonable
            # upper bound to this) iterations have been reached
            run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability_at_this_temperature, maximum_iterations=Int(ceil(relaxation_iterations_vector[temperature_index])), verbose=false, configuration_correlation_convergence_criteria=exp(-1))

            E = energy(cube)
            E_running_total += E
            E_squared_running_total += E^2

            if creating_energy_histogram && insorted(T, energy_histogram_sample_temperatures)
                energy_samples_by_temperature[indexin(T,energy_histogram_sample_temperatures), sample_index] .= E
            end

            if collecting_neighbour_energy_info && insorted(T, neighbour_sample_temperatures)
                # If neighbour_per_configuration_sample_size=0 (default) then just collect all neighbours, otherwise collect a random sample of neighbours
                if neighbours_per_configuration_sample_size==0
                    samples = all_neighbour_energies(cube, collecting_swap_move_neighbours, keep_energy_deltas_only=!collect_initial_and_final_energies, neighbour_moves_away=neighbour_moves_away)
                else
                    samples = sample_neighbour_energies(cube, collecting_swap_move_neighbours, neighbours_per_configuration_sample_size; keep_energy_deltas_only=!collect_initial_and_final_energies, neighbour_moves_away=neighbour_moves_away)
                end

                # If collecting minimum neighbour energy delta only then just store minimum
                if collect_minimum_neighbour_energy_info_only
                    neighbour_information_by_temperature[indexin(T,neighbour_sample_temperatures), (sample_index-1)+1:sample_index] .=  minimum(samples) 
                else
                    neighbour_information_by_temperature[indexin(T,neighbour_sample_temperatures), (sample_index-1)*number_of_neighbours+1:sample_index*number_of_neighbours] .= samples
                end
            end
 
        end

        E_average_by_temperature[temperature_index] = E_running_total/average_sample_size_per_temperature
        E_squared_average_by_temperature[temperature_index] = E_squared_running_total/average_sample_size_per_temperature

        measured_relaxation_iterations_by_temperature[temperature_index] = final_iteration_number/2 # (divide by 2 as relaxation stage uses 2 relaxation_iterations before stopping)
        accepted_candidates_by_temperature[temperature_index] = final_accepted_candidates_number/2 # (divide by 2 as relaxation stage uses 2 relaxation_iterations before stopping)
        final_configuration_correlation_function_by_temperature[temperature_index] = final_configuration_correlation_function # This should be e^-2 for relaxed anneal with out current set-up

        # Print current temperature and average energy if verbose_annealing mode activated
        if verbose_annealing
            printstyled("Currently at Temperature:  $T [$(temperature_index)/$(length(temperature_vector))] (P_swap=$swap_move_probability, T_swap=$T_swap, L=$(cube.L))\n"; underline=true)

            if T <= T_swap && swap_move_probability!=0.0
                println("Using swap moves at this temperature")
            end

            println("Average Energy: $(E_average_by_temperature[temperature_index])")
            println("-Average Energy/Solved Configuration Energy: $(-E_average_by_temperature[temperature_index]/solved_configuration_energy(cube))")
            println("Relaxation Iterations to e^-2 Correlation Function: $final_iteration_number")
            println("Accepted Candidates: $final_accepted_candidates_number")
            println("Acceptance Rate: $((final_accepted_candidates_number/final_iteration_number)*100) %")
            println("Relaxation Converged?: $relaxation_converged")
            println("Final Configuration Correlation Function (for t=2*tau): $final_configuration_correlation_function")
        end

        # Break out of loop if user types 'break
        # Check for user input
        catch e 
            println("Stopping the loop due to user input!")

            # Trim all arrays so curernt temperature is final temperature
            temperature_vector = temperature_vector[1:temperature_index]
            E_average_by_temperature = E_average_by_temperature[1:temperature_index]
            E_squared_average_by_temperature = E_squared_average_by_temperature[1:temperature_index]
            measured_relaxation_iterations_by_temperature = measured_relaxation_iterations_by_temperature[1:temperature_index]
            accepted_candidates_by_temperature = accepted_candidates_by_temperature[1:temperature_index]
            final_configuration_correlation_function_by_temperature = final_configuration_correlation_function_by_temperature[1:temperature_index]

            if collecting_neighbour_energy_info
                neighbour_information_by_temperature = neighbour_information_by_temperature[1:temperature_index,:]
            end

            if creating_energy_histogram
                energy_samples_by_temperature = energy_samples_by_temperature[1:temperature_index,:]
            end

            break
        end
    
    end

    # Return results as dictionary
    if collecting_neighbour_energy_info
        return temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, measured_relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature, neighbour_information_by_temperature
    elseif creating_energy_histogram
        return temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, measured_relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature, energy_samples_by_temperature
    else
        return temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, measured_relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature 
    end
end