include("relaxed_anneal_experiment.jl")


function infinite_temperature_asymptotics(N=1e8; L::Int64=11)
    N = Int(N)
    N_i = 10000

    cube = RubiksCube(L)
    running_total_energy = 0
    running_total_energy_squared = 0
    running_total_configuration_autocorrelation_function = 0
    running_total_configuration_autocorrelation_function_squared = 0
    running_total_slice_rotations = 0

    # Initial randomisation
    for i in 1:N_i
        random_rotate!(cube)
    end

    initial_cube_configuration = deepcopy(cube.configuration)

    for i in 1:N
        random_rotate!(cube)
        energy_sample = energy(cube)
        running_total_energy += energy_sample
        running_total_energy_squared += energy_sample^2

        configuration_autocorrelation_function_sample = configuration_correlation_function(cube.configuration,initial_cube_configuration)
        running_total_configuration_autocorrelation_function += configuration_autocorrelation_function_sample
        running_total_configuration_autocorrelation_function_squared += configuration_autocorrelation_function_sample^2
        running_total_slice_rotations += 1
    end

    println("Running Total Slice Rotations: ", running_total_slice_rotations)
    println("Final Average Energy Density: ", running_total_energy/(running_total_slice_rotations*solved_configuration_energy(cube)))
    println("Final Average Energy Density Squared: ", running_total_energy_squared/(running_total_slice_rotations*solved_configuration_energy(cube)^2))
    println("Final Average Configuration Autocorrelation Function: ", running_total_configuration_autocorrelation_function/running_total_slice_rotations)
    println("Final Average Configuration Autocorrelation Function Squared: ", running_total_configuration_autocorrelation_function_squared/running_total_slice_rotations)

    average_energy_standard_deviation = sqrt(running_total_energy_squared/(running_total_slice_rotations*solved_configuration_energy(cube)^2) - (running_total_energy/(running_total_slice_rotations*solved_configuration_energy(cube)))^2)
    average_configuration_autocorrelation_function_standard_deviation = sqrt(running_total_configuration_autocorrelation_function_squared/running_total_slice_rotations - (running_total_configuration_autocorrelation_function/running_total_slice_rotations)^2)

    println("Average Energy Standard Deviation: ", average_energy_standard_deviation)
    println("Average Configuration Autocorrelation Function Standard Deviation: ", average_configuration_autocorrelation_function_standard_deviation)

    average_energy_standard_error = average_energy_standard_deviation/sqrt(running_total_slice_rotations)
    average_configuration_autocorrelation_function_standard_error = average_configuration_autocorrelation_function_standard_deviation/sqrt(running_total_slice_rotations)

    println("Average Energy Standard Error: ", average_energy_standard_error)
    println("Average Configuration Autocorrelation Function Standard Error: ", average_configuration_autocorrelation_function_standard_error)

    return running_total_energy, running_total_slice_rotations, cube
end

function large_randomised_cube_autocorrelation_function_asymptotics(L::Int64;N=1e8)

    N = Int(N)
    N_i = 10000

    cube = RubiksCube(L)
    running_total_configuration_autocorrelation_function = 0
    running_total_configuartion_autocorrelation_function_squared = 0
    running_total_slice_rotations = 0

    initial_cube_configuration = deepcopy(cube.configuration)

    # Make randomised cube
    facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
            shuffle!(facelets)
            new_faces = reshape(facelets, 6, L, L)
            for i in 1:6
                initial_cube_configuration[i][:,:] .= new_faces[i,:,:]
    end

    cube.configuration .= initial_cube_configuration
    initial_cube_configuration = deepcopy(cube.configuration)

    for i in 1:N
        random_rotate!(cube)
        configuration_autocorrelation_function_sample = configuration_correlation_function(cube.configuration,initial_cube_configuration)
        running_total_configuration_autocorrelation_function += configuration_autocorrelation_function_sample
        running_total_configuartion_autocorrelation_function_squared += configuration_autocorrelation_function_sample^2

        running_total_slice_rotations += 1
    end

    println("Running Total Slice Rotations: ", running_total_slice_rotations)
    println("Final Average Configuration Autocorrelation Function: ", running_total_configuration_autocorrelation_function/running_total_slice_rotations)
    println("Final Average Configuration Autocorrelation Function Squared: ", running_total_configuartion_autocorrelation_function_squared/running_total_slice_rotations)

    average_configuration_autocorrelation_function_standard_deviation = sqrt(running_total_configuartion_autocorrelation_function_squared/running_total_slice_rotations - (running_total_configuration_autocorrelation_function/running_total_slice_rotations)^2)

    println("Average Configuration Autocorrelation Function Standard Deviation: ", average_configuration_autocorrelation_function_standard_deviation)

    average_configuration_autocorrelation_function_standard_error = average_configuration_autocorrelation_function_standard_deviation/sqrt(running_total_slice_rotations)

    println("Average Configuration Autocorrelation Function Standard Error: ", average_configuration_autocorrelation_function_standard_error)

    return running_total_configuration_autocorrelation_function/running_total_slice_rotations

end



large_randomised_cube_autocorrelation_function_asymptotics(100, N=1e8)

# infinite_temperature_asymptotics(1e9, L=11)