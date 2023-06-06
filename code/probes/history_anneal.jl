include("../core/rubiks_cube.jl")
include("../core/monte_carlo.jl")


function history_anneal!(cube::RubiksCube, temperature_vector::Vector{Float64}; swap_move_probability::Float64=0.0, T_swap::Float64=0.0, verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false, mix::Bool=true)

    # Runs history anneal on Rubik's Cube using the Metroolis+Swap algorithm at each temperature to probe properties of the cube.
    # Unlike relaxed anneal we do not relax the cube at each temperature, and we only take one energy sample per temperature
    # But typically we use a much longer temperature vector to still cool slowly

    # Mixing Stage ---

    if mix

        # (Randomise Rubik's cube initially by running Metropolis algorithm at infinite temperature)
        # i.e. beta = 0, swap_move_probability = 0, and by default have maximum_iterations = 10,000 or a configuration correlation function value of e^-10
        run_metropolis_swap_algorithm!(cube, 0.0; swap_move_probability=0.0, maximum_iterations=10000, verbose=false, configuration_correlation_convergence_criteria=exp(-10))

        if verbose_annealing
            printstyled("New Cube With P_swap = $swap_move_probability below T_swap = $T_swap \n"; color=:blue)
            println("Mixed cube")
            println("Cube Energy/Infinite Temperature Energy: $(energy(cube)/infinite_temperature_energy(cube))")
        end

    end


    # Annealing Stage ----

    # Create arrays to store parameters for each temeprature
    E_by_temperature = zeros(length(temperature_vector))
    
    # Cool Rubik's cube from T_1 to T_0 by temperatures described in the temperature vector
    for (temperature_index, T) in pairs(temperature_vector)
        beta = 1/T

        # Only allow swap moves below a certain temperature T_swap 
        if T <= T_swap && swap_move_probability!=0.0
            swap_move_probability_at_this_temperature = swap_move_probability
        else
            swap_move_probability_at_this_temperature = 0.0
        end

        # Measurement Stage ---

        # Measure E ONCE at this temperature after a single Metropolis+Swap algorithm step
        run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability_at_this_temperature, maximum_iterations=1, verbose=false)


        E_by_temperature[temperature_index] = energy(cube)

        # Print normalised energy values if verbose_annealing mode activated
        if verbose_annealing && mod(temperature_index,(length(temperature_vector)-1)/10) == 0
            printstyled("Currently at Temperature:  $T [$(temperature_index)/$(length(temperature_vector))] (P_swap=$swap_move_probability, T_swap=$T_swap, L=$(cube.L))\n"; underline=true)
            println("-Energy/Solved Configuration Energy: $(-E_by_temperature[temperature_index]/solved_configuration_energy(cube))")
        end
    
    end

    # Return results as dictionary
    return temperature_vector, E_by_temperature

end