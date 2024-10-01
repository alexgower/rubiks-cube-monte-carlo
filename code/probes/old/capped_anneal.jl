include("../core/rubiks_cube.jl")
include("../core/monte_carlo.jl")


function capped_anneal!(cube::RubiksCube, temperature_vector::Vector{Float64}, energy_cap::Float64, energy_floor::Float64; swap_move_probability::Float64=0.0, verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false)

    # Runs capped anneal on Rubik's Cube using the Metroolis+Swap algorithm at each temperature to probe properties of the cube.
    # Unlike relaxed anneal we do not relax the cube at each temperature, and we only take one energy sample per temperature
    # But typically we use a much longer temperature vector to still cool slowly
    # Also we only allow configurations below a certain energy cap to be accepted
    # Also we save all configurations with the energy of the energy floor, and their relative configuration correlation functions to the original configuration


    reference_cube = RubiksCube(cube.L)
    reference_cube.configuration = deepcopy(cube.configuration)

    # Annealing Stage ----

    # Create arrays to store parameters for each temperature step
    E_by_temperature_step = zeros(length(temperature_vector))
    configuration_correlation_function_by_temperature_step = zeros(length(temperature_vector))

    energy_floor_configurations = []
    energy_floor_configuration_correlation_function_values = []
    
    # Run MCMC anneal on Cube using temperatures described in the temperature vector
    for (temperature_index, T) in pairs(temperature_vector)
        beta = 1/T

        if swap_move_probability == 0.0
            candidate_generating_function! = random_rotate!
        else
            candidate_generating_function! = random_swap_move!
        end

        monte_carlo_timestep!(cube, candidate_generating_function!, beta, verbose=verbose_metropolis_swap, energy_cap=energy_cap)


        E_by_temperature_step[temperature_index] = energy(cube)
        configuration_correlation_function_by_temperature_step[temperature_index] = configuration_correlation_function(cube, reference_cube)

        if E_by_temperature_step[temperature_index] == energy_floor
            push!(energy_floor_configurations, deepcopy(cube.configuration))
            push!(energy_floor_configuration_correlation_function_values, configuration_correlation_function_by_temperature_step[temperature_index])
        end

        # Print normalised energy values if verbose_annealing mode activated
        if verbose_annealing && mod(temperature_index,(length(temperature_vector)-1)/10) == 0
            printstyled("Currently at Temperature:  $T [$(temperature_index)/$(length(temperature_vector))] (P_swap=$swap_move_probability), L=$(cube.L))\n"; underline=true)
            println("-Energy/Solved Configuration Energy: $(-E_by_temperature_step[temperature_index]/solved_configuration_energy(cube))")
            println("Absolute Energy: $(E_by_temperature_step[temperature_index])")
            println("Configuration Correlation Function: $(configuration_correlation_function_by_temperature_step[temperature_index])")
        end
    
    end

    # Return results as dictionary
    return temperature_vector, E_by_temperature_step, configuration_correlation_function_by_temperature_step, energy_floor_configurations, energy_floor_configuration_correlation_function_values

end