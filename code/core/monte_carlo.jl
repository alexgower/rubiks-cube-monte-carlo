include("rubiks_cube.jl")
include("swap_moves.jl")


function monte_carlo_timestep!(cube::RubiksCube, candidate_generating_function!::Function, beta::Float64; verbose::Bool=false)
    # Perform Metropolis algorithm Monte Carlo Step on cube and candidate generated configurations at given beta
    # i.e. here we either reverse the candidate_generating_function's action on the cube (i.e. reject candidate configuration) or not

    current_energy = energy(cube)

    if verbose
        println("Current Configuration: ")
        println(cube.configuration)
    end

    # Generate candidate configuration and calculate it's energy
    # This function just returns 'candidate_reversing_information' after modifying cube so we know how to reverse modification
    candidate_reversing_information = candidate_generating_function!(cube)

    candidate_energy = energy(cube)

    if verbose
        println("Candidate Configuraiton: ")
        println(cube.configuration)
        println("Candidate Generating Rotation Information: ")
        println(candidate_reversing_information) 
    end

    # Implement the acceptance distribution A(x'|x) = min(1, P(x')/P(x)) by calculating the acceptance
    # ratio alpha = P(x')/P(x), then generating a random number u between 0 and 1, then accepting the candidate
    # if u <= alpha



    # Only print information if verbose mode on
    if verbose
        println("Current Energy: $current_energy")
        println("Candidate Energy: $candidate_energy")
        println("Alpha: $(exp(beta * (current_energy - candidate_energy)))")
    end


    alpha = exp(beta * (current_energy - candidate_energy))

    if rand() <=  alpha # If the acceptance probability alpha is larger than u then we accept the new state i.e do not reverse it

        if verbose
            printstyled("Switched \n"; color=:green)

            if alpha < 1.0
                printstyled("To a higher energy \n"; color=:red)
            end
        end

        # Also return accepted_candidates_increase = 1
        return 1
    else

        # Otherwise we reject the candidate configuration and revert it to the original configuraiton
        candidate_generating_function!(cube;reverse=true,candidate_reversing_information=candidate_reversing_information)

        # Also return accepted_candidates_increase = 0
        return 0
    end
end





function run_metropolis_swap_algorithm!(cube::RubiksCube, beta::Float64; swap_move_probability::Float64=0.0, maximum_iterations::Int64=1000, configuration_correlation_convergence_criteria::Float64=exp(-1), verbose::Bool=false)
    # Runs local Metropolis+Swap algorithm on given Rubik's cube, with composite swap moves also included with a certain
    # likelihood, until either configuration correlation function converges or maximum number of iterations is reached.

    # Make initial count variable and variable to count accepted candidates
    current_iteration = 0
    accepted_candidates = 0

    # Make initial configuration correlation function convergence scheme parameters
    initial_cube = RubiksCube(cube.L)
    initial_cube.configuration = deepcopy(cube.configuration)
    current_configuration_correlation_function_value = configuration_correlation_function(cube,initial_cube)

    while (current_iteration <= maximum_iterations) && (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)

        # Pick a random number between 0 and 1, 
        # If it is smaller than the swap_move_probability then we propose a swap move candidate configuration, else we propose a slice rotated candidate configuration
        if rand() < swap_move_probability
            candidate_generating_function! = random_swap_move!
            
            if verbose
                printstyled("Proposed swap move \n"; color=:blue)
            end
        else
            # (The random_rotate() function within the Rubik's Cube class just implements proposal probability g(x'|x) = 1/|R_{f,l,o}| for
            # whatever sized Rubik's cube we are dealing with)
            candidate_generating_function! = random_rotate!
        end

        if verbose
            printstyled("Current Iteration: $current_iteration \n", underline=true)
        end

        # Now do a Monte Carlo timestep using this candidate configuration at this beta
        accepted_candidates_increase = monte_carlo_timestep!(cube, candidate_generating_function!, beta; verbose=verbose)

        # Update iteration number, accepted_candidates and configuration correlation function and go to next iteration
        current_iteration += 1
        accepted_candidates += accepted_candidates_increase
        current_configuration_correlation_function_value = configuration_correlation_function(cube, initial_cube)
    end

    # If current_iteration is less than maximum_iterations then we know that the convergence criteria was satisfied
    converged = (current_iteration < maximum_iterations)
    if verbose
        if converged
            printstyled("Configuration Correlation Function Converged!", color=:green)
        else
            printstyled("Maximum Iterations Reached!", color=:red)
        end
    end

    # Return some useful information
    return (converged, current_configuration_correlation_function_value, current_iteration, accepted_candidates)
end