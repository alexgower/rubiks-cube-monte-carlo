using Distributed
using SharedArrays
using Plots
using DelimitedFiles
using LaTeXStrings

include("../probes/autocorrelation_anneal.jl")




@inbounds @fastmath function fast_autocorrelation_experiment(simulation_name::String, L::Int64, annealing_swap_move_probability::Float64, autocorrelation_swap_move_probability::Float64, T_1::Float64, T_0::Float64, N_T::Int64, sample_temperature::Float64, relaxation_iterations_per_temperature::Int64, autocorrelation_window_length::Int64; verbose_annealing::Bool=true, verbose_metropolis_swap::Bool=false, verbose_graph_annealing::Bool=false, original_configuration::Union{Vector{Matrix{Int64}},Nothing}=nothing, inherent_disorder::Bool=false, lag_limit::Int64=100)

    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    all_temperatures = unique(sort(vcat(temperature_vector, sample_temperature), rev=true))

    average_sample_size_per_temperature = 1

    ## -- MAKE INITIAL CUBE --
    cube = RubiksCube(L)

    if !isnothing(original_configuration)
        cube.configuration = original_configuration
    
    elseif inherent_disorder
        facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
        shuffle!(facelets)
        new_faces = reshape(facelets, 6, L, L)
        for i in 1:6
            cube.configuration[i][:,:] .= new_faces[i,:,:]
        end
    end

    println("Using initial cube configuration: $(cube.configuration)")

    ## -- RUN ANNEALING --


    # Runs temperature anneal on Rubik's Cube using the Metropolis+Swap algorithm at each temperature to
    # relax and probe properties of the cube.
    # Runs autocorrelation anneal at sample temperature where measure every MC step for autocorrelation_window_length 
    # to deduce autocorrelation time and stretching exponent and get autocorrelation function decay
    reference_cube = RubiksCube(cube.L)


    # Mixing Stage ---

    run_metropolis_swap_algorithm!(cube, 0.0; swap_move_probability=0.0, maximum_iterations=Int(ceil(10*relaxation_iterations_per_temperature)), verbose=false, configuration_correlation_convergence_criteria=exp(-10))

    if verbose_annealing
        printstyled("New Cube With P_swap = $annealing_swap_move_probability \n"; color=:blue)
        println("Mixed cube")
        println("Cube Energy/Infinite Temperature Energy: $(energy(cube)/infinite_temperature_energy(cube))")
    end



    # Annealing Stage ----

    configuration_autocorrelation_function_by_lags = zeros(lag_limit+1)

    # Cool Rubik's cube from T_1 to T_0 by temperatures described in the temperature vector
    for (temperature_index, T) in pairs(all_temperatures)
        if T < sample_temperature
            break
        end

        println("Annealing at temperature T = $T")

        beta = 1/T

        # ANNEAL WITH annealing_swap_move_probability
        run_metropolis_swap_algorithm!(cube, beta, swap_move_probability=annealing_swap_move_probability, maximum_iterations=relaxation_iterations_per_temperature, verbose=false, configuration_correlation_convergence_criteria=exp(-2))


        if T == sample_temperature
                println("Starting Autocorrelation at temperature T = $T")

                # Do average_sample_size_per_temperature runs
                for sample_index in 1:average_sample_size_per_temperature
                    reference_cube.configuration = deepcopy(cube.configuration)

                    # Create arrays to store parameters for each sample
                    configuration_by_sample = [empty(cube.configuration) for i in 1:autocorrelation_window_length]
    

                    # Do autocorrelation_window_length MC steps and save every energy value and configuration
                    for step_index in 1:autocorrelation_window_length
                        # AUTOCORRELATION WITH autocorrelation_swap_move_probability
                        candidate_generating_function! = autocorrelation_swap_move_probability==0.0 ? random_rotate! : random_swap_move!
                        monte_carlo_timestep!(cube, candidate_generating_function!, beta, verbose=false)

                        configuration_by_sample[step_index] = deepcopy(cube.configuration)
                    end

                    lags = collect(0:lag_limit) # i.e. calculate for tau from 0 to 100 MC steps behind

                    configuration_autocorrelation_function_by_lags = configuration_autocorrelation(configuration_by_sample,lags)
                end

            end
    end

    println("Finished Annealing")


    ### --- SAVE RESULTS ---

    ## -- SAVE CONFIGURATION AVERAGES BY TIME --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_configuration_autocorrelation_averages_by_time.csv"

    touch(joinpath("results/autocorrelation_anneal_results/fast_results",filename))

    open(joinpath("results/autocorrelation_anneal_results/fast_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, Samples in Average N, Configuration autocorrelation by lag \n")
    
        # Use join to convert array values to a comma-separated string
        str = join(configuration_autocorrelation_function_by_lags, ", ")
        write(simulation_file, "$(sample_temperature), $(average_sample_size_per_temperature), $str \n")

    end

end

function configuration_autocorrelation(configuration_by_sample, lags::Vector{Int64})
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