using Distributed
using SharedArrays
using Plots
using DelimitedFiles
using LaTeXStrings

include("../probes/autocorrelation_anneal.jl")




@inbounds @fastmath function autocorrelation_experiment(simulation_name::String, L::Int64, annealing_swap_move_probability::Float64, autocorrelation_swap_move_probability::Float64, T_1::Float64, T_0::Float64, N_T::Int64, sample_temperatures::Vector{Float64}, relaxation_iterations_per_temperature::Int64, average_sample_size_per_temperature::Int64, autocorrelation_window_length::Int64; verbose_annealing::Bool=true, verbose_metropolis_swap::Bool=false, verbose_graph_annealing::Bool=false, original_configuration::Union{Vector{Matrix{Int64}},Nothing}=nothing, inherent_disorder_average::Bool=false, parallel_anneals::Int64=1, lag_limit::Int64=100)

    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]

    samples_per_temperature_per_anneal = Int(ceil(average_sample_size_per_temperature/parallel_anneals))


    # Initialize SharedArrays with the correct sizes
    energy_autocorrelation_times_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), average_sample_size_per_temperature))
    energy_stretching_exponents_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), average_sample_size_per_temperature))

    configuration_autocorrelation_times_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), average_sample_size_per_temperature))
    configuration_stretching_exponents_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), average_sample_size_per_temperature))

    # Just bin together the energy autocorrelation averages by time rather than saving all runs separately
    energy_autocorrelation_averages_by_time_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), lag_limit+1))
    configuration_autocorrelation_averages_by_time_by_temperature = SharedArray{Float64,2}((length(sample_temperatures), lag_limit+1))

    @sync @distributed for trial in 1:parallel_anneals

        printstyled("Trial: $trial \n", color=:light_blue)

        ## -- MAKE INITIAL CUBE --
        cube = RubiksCube(L)

        if !isnothing(original_configuration)
            cube.configuration = original_configuration
        
        elseif inherent_disorder_average
            facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
            shuffle!(facelets)
            new_faces = reshape(facelets, 6, L, L)
            for i in 1:6
                cube.configuration[i][:,:] .= new_faces[i,:,:]
            end
        end
    
        println("Using initial cube configuration: $(cube.configuration)")

        ## -- RUN ANNEALING --

        _, trial_energy_autocorrelation_time_by_temperature, trial_energy_stretching_exponent_by_temperature, trial_configuration_autocorrelation_time_by_temperature, trial_configuration_stretching_exponent_by_temperature, trial_energy_autocorrelation_average_by_time_by_temperature, trial_configuration_autocorrelation_average_by_time_by_temperature = autocorrelation_anneal!(cube, temperature_vector, relaxation_iterations_per_temperature, sample_temperatures, samples_per_temperature_per_anneal, autocorrelation_window_length, annealing_swap_move_probability=annealing_swap_move_probability, autocorrelation_swap_move_probability=autocorrelation_swap_move_probability, verbose_annealing=verbose_annealing, verbose_metropolis_swap=verbose_metropolis_swap, verbose_graph_annealing=verbose_graph_annealing, lag_limit=lag_limit)


        ## -- STORE IN SHARED ARRAYS --
        
        # Remember all trials have all sample temperatures, and just some of the samples per temperature

        energy_autocorrelation_times_by_temperature[:,(trial-1)*samples_per_temperature_per_anneal+1:trial*samples_per_temperature_per_anneal] .= trial_energy_autocorrelation_time_by_temperature
        energy_stretching_exponents_by_temperature[:,(trial-1)*samples_per_temperature_per_anneal+1:trial*samples_per_temperature_per_anneal] .= trial_energy_stretching_exponent_by_temperature

        configuration_autocorrelation_times_by_temperature[:,(trial-1)*samples_per_temperature_per_anneal+1:trial*samples_per_temperature_per_anneal] .= trial_configuration_autocorrelation_time_by_temperature
        configuration_stretching_exponents_by_temperature[:,(trial-1)*samples_per_temperature_per_anneal+1:trial*samples_per_temperature_per_anneal] .= trial_configuration_stretching_exponent_by_temperature

        energy_autocorrelation_averages_by_time_by_temperature[:, :] .+= trial_energy_autocorrelation_average_by_time_by_temperature
        configuration_autocorrelation_averages_by_time_by_temperature[:, :] .+= trial_configuration_autocorrelation_average_by_time_by_temperature

    end

    # Average the autocorrelation averages by time
    energy_autocorrelation_averages_by_time_by_temperature ./= parallel_anneals
    configuration_autocorrelation_averages_by_time_by_temperature ./= parallel_anneals


    println("Finished Annealing")


    ### --- SAVE RESULTS ---

    ### --- SAVE PARAMETERS ---

    ## -- SAVE ENERGY TAUS FILE --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_energy_taus.csv"

    touch(joinpath("results/autocorrelation_anneal_results",filename))

    open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder_average
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, tau_E samples \n")
        
        for temperature_index in 1:length(sample_temperatures)
            # Use join to convert array values to a comma-separated string
            autocorr_times_str = join(energy_autocorrelation_times_by_temperature[temperature_index, :], ", ")
            write(simulation_file, "$(sample_temperatures[temperature_index]), $autocorr_times_str \n")
        end

    end


    ## -- SAVE ENERGY BETAS FILE --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_energy_betas.csv"

    touch(joinpath("results/autocorrelation_anneal_results",filename))

    open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder_average
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, beta_E samples \n")
        
        for temperature_index in 1:length(sample_temperatures)
            # Use join to convert array values to a comma-separated string
            str = join(energy_stretching_exponents_by_temperature[temperature_index, :], ", ")
            write(simulation_file, "$(sample_temperatures[temperature_index]), $str \n")
        end

    end


    ## -- SAVE CONFIGURATION TAUS FILE --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_configuration_taus.csv"

    touch(joinpath("results/autocorrelation_anneal_results",filename))

    open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder_average
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, tau_C samples \n")
        
        for temperature_index in 1:length(sample_temperatures)
            # Use join to convert array values to a comma-separated string
            str = join(configuration_autocorrelation_times_by_temperature[temperature_index, :], ", ")
            write(simulation_file, "$(sample_temperatures[temperature_index]), $str \n")
        end

    end


    ## -- SAVE CONFIGURATION BETAS FILE --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_configuration_betas.csv"

    touch(joinpath("results/autocorrelation_anneal_results",filename))

    open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder_average
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, beta_C samples \n")
        
        for temperature_index in 1:length(sample_temperatures)
            # Use join to convert array values to a comma-separated string
            str = join(configuration_stretching_exponents_by_temperature[temperature_index, :], ", ")
            write(simulation_file, "$(sample_temperatures[temperature_index]), $str \n")
        end

    end





    ### --- SAVE AVERAGES BY TIME ---

    ## -- SAVE ENERGY AVERAGES BY TIME --
    filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_energy_autocorrelation_averages_by_time.csv"

    touch(joinpath("results/autocorrelation_anneal_results",filename))

    open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
        write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
        if !isnothing(original_configuration)
            write(simulation_file, "Original Configuration = $original_configuration \n")
        else
            if inherent_disorder_average
                write(simulation_file, "Inherent Disorder Average \n")
            else
                write(simulation_file, "Solved Initial Configuration \n")
            end
        end
        write(simulation_file, "Sample Temperature T, Samples in Average N, Energy autocorrelation average by lag \n")
        
        for temperature_index in 1:length(sample_temperatures)
            # Use join to convert array values to a comma-separated string
            str = join(energy_autocorrelation_averages_by_time_by_temperature[temperature_index, :], ", ")
            write(simulation_file, "$(sample_temperatures[temperature_index]), $(average_sample_size_per_temperature), $str \n")
        end



        ## -- SAVE CONFIGURATION AVERAGES BY TIME --
        filename = simulation_name * '_' * string(annealing_swap_move_probability) * "_configuration_autocorrelation_averages_by_time.csv"

        touch(joinpath("results/autocorrelation_anneal_results",filename))

        open(joinpath("results/autocorrelation_anneal_results",filename), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$annealing_swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelation_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
            if !isnothing(original_configuration)
                write(simulation_file, "Original Configuration = $original_configuration \n")
            else
                if inherent_disorder_average
                    write(simulation_file, "Inherent Disorder Average \n")
                else
                    write(simulation_file, "Solved Initial Configuration \n")
                end
            end
            write(simulation_file, "Sample Temperature T, Samples in Average N, Configuration autocorrelation average by lag \n")
            
            for temperature_index in 1:length(sample_temperatures)
                # Use join to convert array values to a comma-separated string
                str = join(configuration_autocorrelation_averages_by_time_by_temperature[temperature_index, :], ", ")
                write(simulation_file, "$(sample_temperatures[temperature_index]), $(average_sample_size_per_temperature), $str \n")
            end

        end

    end

end
