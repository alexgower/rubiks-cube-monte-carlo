using FFTW
using LsqFit
using StatsBase
using LaTeXStrings


include("../core/rubiks_cube.jl")
include("../core/monte_carlo.jl")





@inbounds @fastmath function autocorrelation_anneal!(cube::RubiksCube, temperature_vector::Vector{Float64}, relaxation_iterations_per_temperature::Int64=10000, sample_temperatures::Vector{Float64}=empty([0.0]), average_sample_size_per_temperature::Int64=10, autocorrelation_window_length::Int64=1000; swap_move_probability::Float64=0.0,  verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false)

    # Notes ---

    # Runs temperature anneal on Rubik's Cube using the Metropolis+Swap algorithm at each temperature to
    # relax and probe properties of the cube.
    # Runs autocorrelation anneals at each sample temperature where measure every MC step for autocorrelation_window_length 
    # to deduce autocorrelation time and stretching exponent

    sample_temperatures = sort(sample_temperatures, rev=true)

    all_temperatures = sort(vcat(temperature_vector, sample_temperatures), rev=true)

    reference_cube = RubiksCube(cube.L)


    # Mixing Stage ---

        run_metropolis_swap_algorithm!(cube, 0.0; swap_move_probability=0.0, maximum_iterations=Int(ceil(10*relaxation_iterations_per_temperature)), verbose=false, configuration_correlation_convergence_criteria=exp(-10))

        if verbose_annealing
            printstyled("New Cube With P_swap = $swap_move_probability \n"; color=:blue)
            println("Mixed cube")
            println("Cube Energy/Infinite Temperature Energy: $(energy(cube)/infinite_temperature_energy(cube))")
        end



    # Annealing Stage ----

    # Create arrays to store parameters for each temeprature

    energy_autocorrelation_time_by_temperature = zeros(length(sample_temperatures), average_sample_size_per_temperature)
    energy_stretching_exponent_by_temperature = zeros(length(sample_temperatures), average_sample_size_per_temperature)

    configuration_autocorrelation_time_by_temperature = zeros(length(sample_temperatures), average_sample_size_per_temperature)
    configuration_stretching_exponent_by_temperature = zeros(length(sample_temperatures), average_sample_size_per_temperature)

    
    # Cool Rubik's cube from T_1 to T_0 by temperatures described in the temperature vector
    for (temperature_index, T) in pairs(all_temperatures)
        # try
            beta = 1/T

            run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=swap_move_probability, maximum_iterations=relaxation_iterations_per_temperature, verbose=false, configuration_correlation_convergence_criteria=exp(-2))


            if T ∈ sample_temperatures

                    # Do average_sample_size_per_temperature runs
                    for sample_index in 1:average_sample_size_per_temperature
                        reference_cube.configuration = deepcopy(cube.configuration)

                        # Create arrays to store parameters for each sample
                        energy_by_sample = zeros(autocorrelation_window_length)
                        energy_autocorrelation_function_by_lags = zeros(autocorrelation_window_length)
                        configuration_by_sample = []
        

                        # Do autocorrelation_window_length MC steps and save every energy correlation function value and relative configuration correlation function value 
                        for step_index in 1:autocorrelation_window_length
                            candidate_generating_function! = swap_move_probability==0.0 ? random_rotate! : random_swap_move!
                            monte_carlo_timestep!(cube, candidate_generating_function!, beta, verbose=false)

                            energy_by_sample[step_index] = energy(cube)
                            push!(configuration_by_sample, deepcopy(cube.configuration))
                        end

                        lags = collect(0:100) # i.e. calculate for tau from 0 to 100 MC steps behind

                        energy_autocorrelation_function_by_lags = autocor(energy_by_sample,lags)
                        configuration_autocorrelation_function_by_lags = configuration_autocorrelation(configuration_by_sample,lags)

                        println("Energy Autocorrelation Function: $(energy_autocorrelation_function_by_lags)")
                        println("Configuration Autocorrelation Function: $(configuration_autocorrelation_function_by_lags)")

                        # Plot energy autocorrelation function if verbose_annealing mode activated
                        if verbose_annealing
                            graph = plot(energy_autocorrelation_function_by_lags, xlabel="Time [MC Steps]", ylabel="Energy Autocorrelation Function, "*L"\langle E(0)E(t) \rangle", title="Rubik's Cube Anneal, L=$(cube.L), T=$T", label="Energy Autocorrelation Function, "*L"\langle E(0)E(t) \rangle")
                            plot!(graph, configuration_autocorrelation_function_by_lags, label="Configuration Autocorrelation Function, "*L"\langle C(0)C(t) \rangle")
                            display(graph)
                        end

                        # Fit for autocorrelation time and stretching exponent for energy correlation function
                        energy_autocorrelation_model(t, p) = exp.(-(t./p[1]).^p[2])
                        # Fit for autocorrelation time and stretching exponent for configuration correlation function
                        # Remember it asymptotes at 1/6 not 0 and still starts at 1
                        configuration_autocorrelation_model(t, p) = (5/6).* exp.(-(t./p[1]).^p[2]) .+ (1/6)

                        p0 = [10, 0.9]
                        lb = [1e-3, 0.1] # example lower bounds
                        ub = [1e10, 2.0]  # example upper bounds
                        t_values = lags
                        
                        energy_autocorrelation_fit = curve_fit(energy_autocorrelation_model, t_values, energy_autocorrelation_function_by_lags, p0, lower=lb, upper=ub)
                        configuration_autocorrelation_fit = curve_fit(configuration_autocorrelation_model, t_values, configuration_autocorrelation_function_by_lags, p0, lower=lb, upper=ub)

                        energy_autocorrelation_time_by_temperature[indexin(T,sample_temperatures),sample_index] .= energy_autocorrelation_fit.param[1]
                        energy_stretching_exponent_by_temperature[indexin(T,sample_temperatures),sample_index] .= energy_autocorrelation_fit.param[2]
                        configuration_autocorrelation_time_by_temperature[indexin(T,sample_temperatures),sample_index] .= configuration_autocorrelation_fit.param[1]
                        configuration_stretching_exponent_by_temperature[indexin(T,sample_temperatures),sample_index] .= configuration_autocorrelation_fit.param[2]

                        # Print autocorrelation time and stretching exponent if verbose_annealing mode activated
                        if verbose_annealing
                            printstyled("Currently at Sample Temperature:  $T [$(temperature_index)/$(length(temperature_vector))] (P_swap=$swap_move_probability, L=$(cube.L))\n"; underline=true)
                            println("Energy Autocorrelation Time: $(energy_autocorrelation_fit.param[1])")
                            println("Energy Stretching Exponent: $(energy_autocorrelation_fit.param[2])")
                            println("Configuration Autocorrelation Time: $(configuration_autocorrelation_fit.param[1])")
                            println("Configuration Stretching Exponent: $(configuration_autocorrelation_fit.param[2])")
                        end
    
                    end
    
                end

        # Break out of loop if user interrupts
        # Check for user input
        # catch e 
        #     println("Stopping the loop due to user input!")

        #     # Trim all arrays so curernt temperature is final temperature
        #     temperature_vector = temperature_vector[1:temperature_index]
        #     energy_autocorrelation_time_by_temperature = energy_autocorrelation_time_by_temperature[1:temperature_index,:]
        #     energy_stretching_exponent_by_temperature = energy_stretching_exponent_by_temperature[1:temperature_index,:]
        #     configuration_autocorrelation_time_by_temperature = configuration_autocorrelation_time_by_temperature[1:temperature_index,:]
        #     configuration_stretching_exponent_by_temperature = configuration_stretching_exponent_by_temperature[1:temperature_index,:]

        #     break
        # end
    
    end

    # Return results as dictionary
    return sample_temperatures, energy_autocorrelation_time_by_temperature, energy_stretching_exponent_by_temperature, configuration_autocorrelation_time_by_temperature, configuration_stretching_exponent_by_temperature
end


function configuration_autocorrelation(configuration_by_sample::Vector{Any}, lags::Vector{Int64})
    # Calculate autocorrelation function for configuration
    # configuration_by_sample is a vector of configurations
    # lags is a vector of lags to calculate autocorrelation function for

    # Calculate autocorrelation function for configuration
    configuration_autocorrelation_function = zeros(length(lags))
    for (lag_index, lag) in pairs(lags)
        configuration_autocorrelation_function[lag_index] = mean([configuration_correlation_function(configuration_by_sample[i], configuration_by_sample[i+lag]) for i in 1:length(configuration_by_sample)-lag])
    end

    return (1/configuration_autocorrelation_function[1]) .* configuration_autocorrelation_function
end