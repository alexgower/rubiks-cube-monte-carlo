# ----------

using Combinatorics

function test_move_reversals(cube::RubiksCube)

    number_of_three_cycles = 0
    number_of_opposite_orientation_rotations = 0

    println("ORIENTATION SWAP MOVES TEST ---")

    if isodd(cube.L)
        cubelet_subsystem_labels = ["sigma_", "tau_"]
    else
        cubelet_subsystem_labels = ["sigma_"]
    end

    for cubelet_subsystem_label in cubelet_subsystem_labels

        if cubelet_subsystem_label == "sigma_"
            number_of_cubelets_in_subsystem = 8
        else # (cubelet_subsystem_label=="tau_" case)
            number_of_cubelets_in_subsystem = 12
        end

        all_random_cubelet_indices_combinations = permutations(1:number_of_cubelets_in_subsystem, 2)

        number_of_opposite_orientation_rotations += length(all_random_cubelet_indices_combinations)

        for random_cubelet_indices_combination in all_random_cubelet_indices_combinations
            # Note this only tests 1 unit of rotation (as that is all is needed to check working correctly) eventhough corners can take 2 units of rotation too
            

            # println(cube.configuration)

            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2])

            # println(cube.configuration)

            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[2], random_cubelet_indices_combination[1])

            # println(cube.configuration)
            
            if cube.configuration != solved_configuration(cube.L)
                println("FAILURE!")
            end

        end
    end

    println("3-CYCLE SWAP MOVES TEST ---")

    for cubelet_subsystem_label in cube.cubelet_subsystems_labels

        if cubelet_subsystem_label=="sigma_"
            number_of_cubelets_in_subsystem = 8
        elseif cubelet_subsystem_label=="tau_"
            number_of_cubelets_in_subsystem = 12
        else
            number_of_cubelets_in_subsystem = 24
        end

        all_random_cubelet_indices_combinations = permutations(1:number_of_cubelets_in_subsystem, 3)

        number_of_three_cycles += length(all_random_cubelet_indices_combinations)


        for random_cubelet_indices_combination in all_random_cubelet_indices_combinations

            # println(cube.configuration)

            three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2], random_cubelet_indices_combination[3])

            # println(cube.configuration)

            three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[3], random_cubelet_indices_combination[2])

            # println(cube.configuration)
            
            if cube.configuration != solved_configuration(cube.L)
                println("FAILURE!")
            end

        end
    end


    println("TESTS COMPLETED SUCCESSFULLY ---")
    println("Number of Three Cycles: $number_of_three_cycles")
    println("Number of Opposite Orientation Rotations: $number_of_opposite_orientation_rotations")

end

# ----------




