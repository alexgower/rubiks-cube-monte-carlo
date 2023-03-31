include("rubiks_cube.jl")

function anneal!(cube::RubiksCube, temperature_vector::Vector{Float64}; swap_move_probability::Float64=0.0, T_swap::Float64=0.0, relaxation_iterations_vector=nothing, average_sample_size::Int64=100, verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false, relaxation_iterations_finder_mode::Bool=false)

    # Runs temperature anneal on Rubik's Cube using the Metropolis+Swap algorithm at each temperature to
    # relax and probe properties of the cube.

    # Notes:
    # - relaxation_iterations_vector: vector of maximum_iterations to use for each temperature in the temperature vector
    # [default will be tau(T) = tau(T_1)*[(tau(T_0)/tau(T_1)]^((T1-T)/(T_1-T_0))]
    # [where tau(T_1) = |R_{f,l,o}| = number of generators of Rubik's rotation group and tau(T_0) = big_tau]
    # - average_sample_size: number of E and E^2 samples to taken in average calculations for each temperature.
    # - relaxation_iterations_finder_mode: If true, will use flat maximum iterations of big_tau for every temperature, and will record
    # and return the number of iterations needed to relax at each temperature


    # Validation --- 

    # Make sure relaxation_iterations_vector (if it exists) has same number of elements as temperature_vector
    if !isnothing(relaxation_iterations_vector) && length(relaxation_iterations_vector) != length(temperature_vector)
        throw(ArgumentError("relaxation_iterations_vector must be same length as temperature_vector"))
    end

    # Extract T_1 and T_0 for future uses
    T1 = temperature_vector[1]
    T0 = temperature_vector[end]

    big_tau = 5000

    # Create relaxation_iterations_vector
    if isnothing(relaxation_iterations_vector)
        if relaxation_iterations_finder_mode==false 

            # Set default relaxation_iterations_vector if none provided and relaxation_iterations_finder_mode not on
            tau_0 = big_tau

            # Number of generators of size-L Rubik's cube is number of faces * number of layers per face * number of
            # rotation orientations per layer = 6 * ceil((L-1)/2) * 2
            tau_1 = 6 * ceil((cube.L-1)/2) * 2

            relaxation_iterations_vector = [tau_1 * (tau_0/tau_1)^((T1-T)/(T1-T0)) for T in temperature_vector]

        elseif relaxation_iterations_finder_mode==true

            # If relaxation_iterations_finder_mode is on, set max_iterations as 100,000 as a constant
            relaxation_iterations_vector = [big_tau for T in temperature_vector]
            tau_0 = big_tau
            tau_1 = big_tau
        end

    else    
        tau_0 = relaxation_iterations_vector[1]
        tau_1 = relaxation_iterations_vector[end]

    end



    # Mixing Stage ---

    # (Randomise Rubik's cube initially by running Metropolis+Swap algorithm at infinite temperature)
    # i.e. beta = 0, swap_move_probability = 0, and by default have maximum_iterations = 10*tau_1
    # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with t=0
    # configuration) has dropped to e^(-10) (i.e. 10 relaxation times) or 10*tau_1 (which should be a reasonable upper bound to this) 
    # iterations have been reached
    # TODO restore no swap moves
    run_metropolis_swap_algorithm!(cube, 0.0; swap_move_probability=0.0, maximum_iterations=10*tau_1, verbose=false, configuration_correlation_convergence_criteria=exp(-10))

    # if verbose_annealing
    if true # TODO restore
        println("Mixed cube")
        println("Cube Energy/Infinite Temperature Energy: $(energy(cube)/infinite_temperature_energy(cube))")
        println(cube.configuration)
    end




    # Annealing Stage ----

    # Create arrays to store parameters for each temeprature
    E_average_by_temperature = zeros(length(temperature_vector))
    E_squared_average_by_temperature = zeros(length(temperature_vector))
    relaxation_iterations_by_temperature = zeros(length(temperature_vector))
    accepted_candidates_by_temperature = zeros(length(temperature_vector))
    final_configuration_correlation_function_by_temperature = zeros(length(temperature_vector))

    
    # Cool Rubik's cube from T_1 to T_0 by temperatures described in the temperature vector
    for (temperature_index, T) in pairs(temperature_vector)
        beta = 1/T

        if verbose_annealing
            printstyled("Currently at Temperature:  $T [$(temperature_index)/$(length(temperature_vector))] (P_swap=$swap_move_probability, T_swap=$T_swap, L=$(cube.L))\n"; underline=true)
        end




        # Relaxation Stage ---

        # Run Metropolis+Swap algorithm to relax Rubik's cube at this temperature
        # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with t=0
        # configuration) has dropped to e^(-2) (i.e. 2 relaxation times) or 2*relation_iterations have been reached

        # Only allow swap moves below a certain temperature T_swap 
        if T <= T_swap && swap_move_probability!=0.0
            if verbose_annealing
                println("Using swap moves at this temperature")
            end
            swap_move_probability_at_this_temperature = swap_move_probability
        else
            swap_move_probability_at_this_temperature = 0.0
        end

        relaxation_converged, final_configuration_correlation_function, final_iteration_number, final_accepted_candidates_number = run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability_at_this_temperature, maximum_iterations=2*relaxation_iterations_vector[temperature_index], verbose=verbose_metropolis_swap, configuration_correlation_convergence_criteria=exp(-2))




        # Measurement Stage ---

        # Calculate <E> and <E^2> at this temperature but only take measurements every relation_iterations steps to ensure statistical independence
        E_running_total = 0.0
        E_squared_running_total = 0.0

        for sample_index in 1:average_sample_size

            # Metropolis+Swap algorithm will terminate when either the configuration correlation function (compared with
            # t=0 configuration) has dropped to e^(-1) (i.e. 1 relaxation time) or tau(T) (which should be a reasonable
            # upper bound to this) iterations have been reached
            run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability_at_this_temperature, maximum_iterations=relaxation_iterations_vector[temperature_index], verbose=false, configuration_correlation_convergence_criteria=exp(-1))

            E_running_total += energy(cube)
            E_squared_running_total += energy(cube)^2
        end

        E_average_by_temperature[temperature_index] = E_running_total/average_sample_size
        E_squared_average_by_temperature[temperature_index] = E_squared_running_total/average_sample_size

        if relaxation_iterations_finder_mode
            relaxation_iterations_by_temperature[temperature_index] = final_iteration_number/2 # (divide by 2 as relaxation stage uses 2 relaxation_iterations before stopping)
            accepted_candidates_by_temperature[temperature_index] = final_accepted_candidates_number/2 # (divide by 2 as relaxation stage uses 2 relaxation_iterations before stopping)
            final_configuration_correlation_function_by_temperature[temperature_index] = final_configuration_correlation_function # (these should be around e^-2 not e^-1) TODO think
        end

        # Print current temperature and average energy if verbose_annealing mode activated
        if verbose_annealing
            println("Average Energy: $(E_average_by_temperature[temperature_index])")
            println("-Average Energy/Solved Configuration Energy: $(-E_average_by_temperature[temperature_index]/solved_configuration_energy(cube))")
            println("2*Relaxation Iterations: $final_iteration_number")
            println("Accepted Candidates: $final_accepted_candidates_number")
            println("Acceptance Rate: $((final_accepted_candidates_number/final_iteration_number)*100) %")
            println("Relaxation Converged?: $relaxation_converged")
            println("Final Configuration Correlation Function (for t=2*tau): $final_configuration_correlation_function")
        end
    
    end

    # Return results as dictionary
    return temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature 

end
