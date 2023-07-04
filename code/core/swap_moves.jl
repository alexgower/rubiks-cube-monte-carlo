# Note that we do have 6/4 necessary and sufficient constraints for the configurations of an odd/even Rubik's Cube to be valid (solvable).
# These couple together the permutation group signatures of the various cubelet groups, and make constraints on the cubelet orientations.
# In theory we could use these constraints to validate any Rubik's Cube candidate configuration, but this would be allowing swap moves
# to ALL other Rubik's Cube configurations from the current configuration i.e. we would lose any notion of locality, and also it would
# likely be incredibly inefficient as most candidate configurations would be completely random and likely rejected

# Instead we just use the fact that every cubelet subsystem commuting subgroup C_X of G is isomorphic to an alternating groups A_n 
# (where n is the number of cubelets in subsystem X)
# (i.e. G is the group of all Rubik's Cube rotations) 
# (i.e. C_X is the group of composite rotations which alter the positions of cubelets within a single subsystem X but leave the rest of the cube unchanged)
# This means that they cannot change the signature of the cubelet subsystem permutation
# ***This means that the minimum allowed cubelet configuration permutation change is a 3-cycle within a cubelet subsystem****

# We also use the constraint that the sum of orientations of cubelets in a cubelet subsystem X is fixed
# ***This means that the minimum allowed cubelet configuration orientation change is opposite rotations of 2 cubelet's orientations***

using StatsBase
using Combinatorics

include("rubiks_cube.jl")

@inline function get_cubelet_subsystem(cube::RubiksCube,cubelet_subsystem_label::String)
    # Function to return an array of references to cube's configuration for the cubelets in subsystem_label

    # First chop cubelet_subsystem_label X_i,j into subsystem_type = X and subsystem_details= [i,j]
    subsystem_type = split(cubelet_subsystem_label,"_")[1]
    subsystem_details = split(split(cubelet_subsystem_label,"_")[2], ",")

    # Validation ---

    # (only odd L cubes have central edge or +-centre cubelets)
    if (subsystem_type=="tau" && iseven(cube.L)) || (subsystem_type=="eta" && iseven(cube.L))
        throw(ArgumentError("Only odd L cubes have central edge or +-centre cubelets "))
    end

    # Calculation ---

    subsystem_array = []

    if subsystem_type == "sigma"
        # Note each corner cubelet contains 3 facelets so each element in the subsystem_array will be a 3x1 array of 3 facelet views
        # The order of facelets in this 3x1 array will be in an ANTICLOCKWISE order (this is well defined if you look at a Rubik's Cube)
        # The precise definitions are shown in the supplementary document
        # Note we require anticlockwise facelet order here so can choose to rotate clockwise or anticlockwise later
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        subsystem_array = [[view(cube.configuration[1],1,1), view(cube.configuration[2],cube.L,1), view(cube.configuration[5],1,cube.L)],
                            [view(cube.configuration[1],1,cube.L), view(cube.configuration[6],1,1), view(cube.configuration[2],cube.L,cube.L)],
                            [view(cube.configuration[1],cube.L,1), view(cube.configuration[5],cube.L,cube.L), view(cube.configuration[4],1,1)],
                            [view(cube.configuration[1],cube.L,cube.L), view(cube.configuration[4],1,cube.L), view(cube.configuration[6],cube.L,1)],
                            [view(cube.configuration[3],cube.L,1), view(cube.configuration[5],1,1), view(cube.configuration[2],1,1)],
                            [view(cube.configuration[3],cube.L,cube.L), view(cube.configuration[2],1,cube.L), view(cube.configuration[6],1,cube.L)],
                            [view(cube.configuration[3],1,1), view(cube.configuration[4],cube.L,1), view(cube.configuration[5],cube.L,1)],
                            [view(cube.configuration[3],1,cube.L), view(cube.configuration[6],cube.L,cube.L),view(cube.configuration[4],cube.L,cube.L)]]
        
    elseif subsystem_type == "theta"
        # Note each wing edge cubelet contains 2 facelets so each element in the subsystem_array will be a 2x1 array of 2 facelet views
        # The order of facelets in this 2x1 array will be in an aribtrary but well defined order (as shown in supplementary document)
        # Note we don't require a specific facelet order here as will only need to flip order later
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        i = parse(Int64,subsystem_details[1])

        subsystem_array = [[view(cube.configuration[3],cube.L,1+i), view(cube.configuration[2],1,1+i)],
                            [view(cube.configuration[2],1,cube.L-i),view(cube.configuration[3],cube.L,cube.L-i)],
                            [view(cube.configuration[2],cube.L,1+i), view(cube.configuration[1],1,1+i)],
                            [view(cube.configuration[1],1,cube.L-i), view(cube.configuration[2],cube.L,cube.L-i)],
                            [view(cube.configuration[1],cube.L,1+i), view(cube.configuration[4],1,1+i)],
                            [view(cube.configuration[4],1,cube.L-i), view(cube.configuration[1],cube.L,cube.L-i)],
                            [view(cube.configuration[4],cube.L,1+i), view(cube.configuration[3],1,1+i)],
                            [view(cube.configuration[3],1,cube.L-i), view(cube.configuration[4],cube.L,cube.L-i)],
                            [view(cube.configuration[1],1+i,1), view(cube.configuration[5],1+i,cube.L)],
                            [view(cube.configuration[5],cube.L-i,cube.L), view(cube.configuration[1],cube.L-i,1)],
                            [view(cube.configuration[6],1+i,1), view(cube.configuration[1],1+i,cube.L)],
                            [view(cube.configuration[1],cube.L-i,cube.L), view(cube.configuration[6],cube.L-i,1)],
                            [view(cube.configuration[2],1+i,1), view(cube.configuration[5],1,1+i)],
                            [view(cube.configuration[5],1,cube.L-i), view(cube.configuration[2],cube.L-i,1)],
                            [view(cube.configuration[2],cube.L-i,cube.L), view(cube.configuration[6],1,1+i)],
                            [view(cube.configuration[6],1,cube.L-i), view(cube.configuration[2],1+i,cube.L)],
                            [view(cube.configuration[6],cube.L,1+i), view(cube.configuration[4],1+i,cube.L), ],
                            [view(cube.configuration[4],cube.L-i,cube.L), view(cube.configuration[6],cube.L,cube.L-i)],
                            [view(cube.configuration[4],1+i,1), view(cube.configuration[5],cube.L,cube.L-i)],
                            [view(cube.configuration[5],cube.L,1+i), view(cube.configuration[4],cube.L-i,1)],
                            [view(cube.configuration[3],cube.L-i,cube.L), view(cube.configuration[6],1+i,cube.L)],
                            [view(cube.configuration[6],cube.L-i,cube.L), view(cube.configuration[3],1+i,cube.L)],
                            [view(cube.configuration[5],1+i,1), view(cube.configuration[3],cube.L-i,1)],
                            [view(cube.configuration[3],1+i,1), view(cube.configuration[5],cube.L-i,1)]]
    
    elseif subsystem_type == "tau"
        # Note each central edge cubelet contains 2 facelets so each element in the subsystem_array will be a 2x1 array of 2 facelet views
        # The order of facelets in this 2x1 array will be in an aribtrary but well defined order (as shown in supplementary document)
        # Note we don't require a specific facelet order here as only flips are possible anyway
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        subsystem_array = [[view(cube.configuration[1],1,Int((cube.L+1)/2)), view(cube.configuration[2],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[6],Int((cube.L+1)/2),1), view(cube.configuration[1],Int((cube.L+1)/2),cube.L)],
                            [view(cube.configuration[1],cube.L,Int((cube.L+1)/2)), view(cube.configuration[4],1,Int((cube.L+1)/2))],
                            [view(cube.configuration[5],Int((cube.L+1)/2),cube.L), view(cube.configuration[1],Int((cube.L+1)/2),1)],
                            [view(cube.configuration[3],cube.L,Int((cube.L+1)/2)), view(cube.configuration[2],1,Int((cube.L+1)/2))],
                            [view(cube.configuration[3],1,Int((cube.L+1)/2)), view(cube.configuration[4],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[2],Int((cube.L+1)/2),cube.L), view(cube.configuration[6],1,Int((cube.L+1)/2))],
                            [view(cube.configuration[6],Int((cube.L+1)/2),cube.L), view(cube.configuration[3],Int((cube.L+1)/2),cube.L)],
                            [view(cube.configuration[2],Int((cube.L+1)/2),1), view(cube.configuration[5],1,Int((cube.L+1)/2))],
                            [view(cube.configuration[5],Int((cube.L+1)/2),1), view(cube.configuration[3],Int((cube.L+1)/2),1)],
                            [view(cube.configuration[4],Int((cube.L+1)/2),1), view(cube.configuration[5],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[4],Int((cube.L+1)/2),cube.L), view(cube.configuration[6],cube.L,Int((cube.L+1)/2))]]    
    
    elseif subsystem_type == "x"
        # Note each X-centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        i = parse(Int64,subsystem_details[1])

        @inbounds @simd for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,1+i)], [view(cube.configuration[face_number],1+i,cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,cube.L-i)], [view(cube.configuration[face_number],cube.L-i,1+i)]])
        end

    elseif subsystem_type == "omega"
        # Note each oblique centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)
    
        i = parse(Int64,subsystem_details[1])
        j = parse(Int64,subsystem_details[2])

        @inbounds @simd for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,1+j)], [view(cube.configuration[face_number],1+j,cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,cube.L-j)], [view(cube.configuration[face_number],cube.L-j,1+i)]])
        end

    elseif subsystem_type == "eta"
        # Note each +-centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        i = parse(Int64,subsystem_details[1])

        @inbounds @simd for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,Int((cube.L+1)/2))], [view(cube.configuration[face_number],Int((cube.L+1)/2),cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,Int((cube.L+1)/2))], [view(cube.configuration[face_number],Int((cube.L+1)/2),1+i)]])
        end

    end
    
    return subsystem_array

end





@inline function three_cycle_cubelets!(cube::RubiksCube, cubelet_subsystem_label::String, cubelet_index_1::Int64, cubelet_index_2::Int64, cubelet_index_3::Int64)    
    # Conducts P3_{X,(a,b,c)} i.e. a 3-cycle in cubelet subsystem X of cubelets in positions -->1-->2-->3-->

    # Validation ---
    # (all cubelet subsystems have size 24 except for the corners (8) and central edges (12))
    if max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 24 || min(cubelet_index_1, cubelet_index_2, cubelet_index_3) < 1 || (cubelet_subsystem_label=="sigma_" && max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 8 ) || (cubelet_subsystem_label=="tau_" && max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 12 )
        throw(ArgumentError("Cubelet indices are outside range of cubelet subsystem"))
    end

    # Operation ---

    cubelet_subsystem = get_cubelet_subsystem(cube, cubelet_subsystem_label)

    # No matter how many facelets each cubelet in this cubelet_subsystem contains (hence the for loop):
    # Move cubelet 2's value to cubelet 3's position 
    # Move cubelet 1's value to cubelet 2's position 
    # Move cubelet 3's value to cubelet 1's position 
    @inbounds @simd for facelet_index in eachindex(cubelet_subsystem[cubelet_index_1])
        temp = copy(cubelet_subsystem[cubelet_index_3][facelet_index])

        cubelet_subsystem[cubelet_index_3][facelet_index] .= cubelet_subsystem[cubelet_index_2][facelet_index]
        cubelet_subsystem[cubelet_index_2][facelet_index] .= cubelet_subsystem[cubelet_index_1][facelet_index]
        cubelet_subsystem[cubelet_index_1][facelet_index] .= temp
    end

end





@inline function opposite_rotate_cubelets!(cube::RubiksCube, cubelet_subsystem_label::String, cubelet_index_1, cubelet_index_2)    
    # Conducts O_{X,(a,b),o} i.e. an orientation rotation in cubelet subsystem X of factor +1 (anticlockwise) to cubelet in position 1 and -1 (clockwise) to cubelet in position 2
    # e.g. [facelet_1, facelet_2, facelet_3] --> [facelet_3, facelet_1, facelet_2] is 1 (anticlockwise) rotation unit (as we used anticlockwise convention in defining facelets in cubelet array)

    # Note wing edge orientations are actually fully constrained by their positions so we only have corner and central edge cases that can have these orientation swap moves
    # Note corners have 3 facelets therefore can be rotated clockwise or anticlockwise, but switching cubelet_index_1 and cubelet_index_2 allows both possibilities
    # Note central edgees only have 2 facelets therefore all 'rotations' are equivalent to just flipping the facelets within the cubelet, and switching cubelet_index_1 and cubelet_index_2 makes no difference

    # Validation ---

    # Only allow corners or central edges as subystems in this function, and only allow central edges for odd L cubes (as don't exist for even L cubes)
    if !(cubelet_subsystem_label=="sigma_" || cubelet_subsystem_label=="tau_") || (cubelet_subsystem_label=="tau_" && iseven(cube.L))
        throw(ArgumentError("This cubelet subsystem label ($cubelet_subsystem_label) doesn't make sense for this cube, or doesn't make sense to have an opposite rotation swap move"))
    end



    # Operation ---
    
    cubelet_subsystem = get_cubelet_subsystem(cube,cubelet_subsystem_label)

    if cubelet_subsystem_label=="sigma_"
        # Anticlockwise for cubelet_index_1: Move facelet 2's value to facelet 3, facelet 1's value to facelet 2, and facelet 3's value to facelet 1  
        # i.e. [1,2,3] --> [3,1,2]
        temp = copy(cubelet_subsystem[cubelet_index_1][3])

        cubelet_subsystem[cubelet_index_1][3] .= cubelet_subsystem[cubelet_index_1][2]
        cubelet_subsystem[cubelet_index_1][2] .= cubelet_subsystem[cubelet_index_1][1]
        cubelet_subsystem[cubelet_index_1][1] .= temp
        
        # Clockwise for cubelet_index_2: Move facelet 1's value to facelet 3, facelet 2's value to facelet 1, and facelet 3's value to facelet 2
        # i.e. [1,2,3] --> [2,3,1]
        temp = copy(cubelet_subsystem[cubelet_index_2][3])

        cubelet_subsystem[cubelet_index_2][3] .= cubelet_subsystem[cubelet_index_2][1]
        cubelet_subsystem[cubelet_index_2][1] .= cubelet_subsystem[cubelet_index_2][2]
        cubelet_subsystem[cubelet_index_2][2] .= temp


    else # (cubelet_subsystem_type=="tau_" case)

        # Simply flip both cubelets as they only have 2 faces
        # i.e. set facelet 1 to facelet 2's value and facelet 2 to facelet 1's value
        temp = copy(cubelet_subsystem[cubelet_index_1][1])

        cubelet_subsystem[cubelet_index_1][1] .= cubelet_subsystem[cubelet_index_1][2]
        cubelet_subsystem[cubelet_index_1][2] .= temp

        temp = copy(cubelet_subsystem[cubelet_index_2][1])

        cubelet_subsystem[cubelet_index_2][1] .= cubelet_subsystem[cubelet_index_2][2]
        cubelet_subsystem[cubelet_index_2][2] .= temp

    end
end





@inline function two_cycle_cubelets!(cube::RubiksCube, cubelet_subsystem_label::String, cubelet_index_1::Int64, cubelet_index_2::Int64)    
    # Conducts P2_{X,(a,b} i.e. a 2-cycle in cubelet subsystem X of cubelets in positions 1<-->2

    # Validation ---
    # (all cubelet subsystems have size 24 except for the corners (8) and central edges (12))
    if max(cubelet_index_1, cubelet_index_2) > 24 || min(cubelet_index_1, cubelet_index_2) < 1 || (cubelet_subsystem_label=="sigma_" && max(cubelet_index_1, cubelet_index_2) > 8 ) || (cubelet_subsystem_label=="tau_" && max(cubelet_index_1, cubelet_index_2) > 12 )
        throw(ArgumentError("Cubelet indices are outside range of cubelet subsystem"))
    end



    # Operation ---

    cubelet_subsystem = get_cubelet_subsystem(cube, cubelet_subsystem_label)

    # No matter how many facelets each cubelet in this cubelet_subsystem contains (hence the for loop):
    # Move cubelet 1's value to cubelet 2's position 
    # Move cubelet 2's value to cubelet 1's position 
    @inbounds @simd for facelet_index in eachindex(cubelet_subsystem[cubelet_index_1])
        temp = copy(cubelet_subsystem[cubelet_index_2][facelet_index])

        cubelet_subsystem[cubelet_index_2][facelet_index] .= cubelet_subsystem[cubelet_index_1][facelet_index]
        cubelet_subsystem[cubelet_index_1][facelet_index] .= temp
    end

end


@inline function get_coupled_subsystems_to_two_cycle(cube::RubiksCube, two_cycle_sigma, two_cycle_theta_ks)

    subsystems_to_two_cycle::Vector{String} = []

    # Independent subsystems
    if two_cycle_sigma
        append!(subsystems_to_two_cycle, ["sigma_"])
    end

    @inbounds for (k,two_cycle_theta_k) in pairs(two_cycle_theta_ks)
        if two_cycle_theta_k
            append!(subsystems_to_two_cycle, ["theta_$k"])
        end
    end

    # Constraint 1
    if two_cycle_sigma
        append!(subsystems_to_two_cycle, vcat(["x_$i" for i in 1:ceil(Int, (cube.L-3)/2)]))

        if isodd(cube.L)
            append!(subsystems_to_two_cycle, ["tau_"])
        end
    end

    # Constraint 2
    if isodd(cube.L)
        @inbounds for (k,two_cycle_theta_k) in pairs(two_cycle_theta_ks)
            # Use XOR to determine whether need to 2_cycle eta_k
            if two_cycle_sigma ⊻ two_cycle_theta_k
                append!(subsystems_to_two_cycle, ["eta_$k"])
            end
        end
    end

    # Constraint 3
    @inbounds for (i,two_cycle_theta_i) in pairs(two_cycle_theta_ks)
        @inbounds for (j,two_cycle_theta_j) in pairs(two_cycle_theta_ks[1:i-1])
                # Use XOR to determine whether need to 2_cycle omega_ij
                if two_cycle_sigma ⊻ two_cycle_theta_i ⊻ two_cycle_theta_j
                    append!(subsystems_to_two_cycle, ["omega_$i,$j"])
                end
        end
    end

    return subsystems_to_two_cycle

end


@inline function random_parity_sector_switch_coupled_two_cycles!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)    
    
    if !reverse
        # Conducts many simultaneous P2_{X,(a,b)} i.e. 2-cycle in cubelet subsystem X of cubelets in positions 1<-->2 
        # But in such a way that we respect the fundamental constraints of the cube (i.e. so it is still solvable afterwards)

        # The required fundamental constraints on cubelet position permutations are:
        # [constraints in brackets are only required for odd L cubes]
        # 1. sgn(sigma) = sgn(x_i) [= sgn(tau)]
        # 2. [sgn(eta_k) = sgn(sigma)sgn(theta_k)]
        # 3. sgn(omega_ij) = sgn(sigma)sgn(theta_i)sgn(theta_j) (for i != j)

        # So first we just choose which of {sigma, {theta_k}} we want to perform a 2-cycle on (and thus change their permutation parity)
        # Then we use the above constraints to determine which of {{x_i},[tau],[{eta_k]}],omega_ij} we also need to perform
        

        # Choose which {sigma, {theta_k}} we want to perform a 2-cycle on
        two_cycle_sigma = rand([true,false])
        two_cycle_theta_ks = [rand([true,false]) for k in 1:ceil(Int, (cube.L-3)/2)]

        # Get full list of all subsystems which must be two cycled in order to respect fundamental constraints of cube given which of above sigma and theta_ks are being two cycled
        subsystems_to_two_cycle::Vector{String} = get_coupled_subsystems_to_two_cycle(cube,two_cycle_sigma,two_cycle_theta_ks)
                
        # Now get 2 random cubelet indices for each cubelet subsystem to be two_cycled and pass these to two_cycle_cubelets!
        candidate_reversing_information = [[subsystem_name,0,0] for subsystem_name in subsystems_to_two_cycle]

        @inbounds @simd for cubelet_subsystem_candidate_reversing_information in candidate_reversing_information

            if cubelet_subsystem_candidate_reversing_information[1]=="sigma_"
                number_of_cubelets_in_subsystem = 8
            elseif cubelet_subsystem_candidate_reversing_information[1]=="tau_"
                number_of_cubelets_in_subsystem = 12
            else
                number_of_cubelets_in_subsystem = 24
            end

            # Get two non-duplicated indices within the number of cubelets in the subsystem
            cubelet_subsystem_candidate_reversing_information[2],cubelet_subsystem_candidate_reversing_information[3] = sample(1:number_of_cubelets_in_subsystem, 2, replace=false)

            # Pass to two_cycle_cubelets
            two_cycle_cubelets!(cube,cubelet_subsystem_candidate_reversing_information[1],cubelet_subsystem_candidate_reversing_information[2],cubelet_subsystem_candidate_reversing_information[3])

        end

        return candidate_reversing_information

    # Reversing case
    else
        # Just undo all 2_cycles by doing another 2_cycle on the same cubelets

        @inbounds @simd for cubelet_subsystem_candidate_reversing_information in candidate_reversing_information
            # Pass to two_cycle_cubelets
            two_cycle_cubelets!(cube,cubelet_subsystem_candidate_reversing_information[1],cubelet_subsystem_candidate_reversing_information[2],cubelet_subsystem_candidate_reversing_information[3])
        end


    end

end





@inline @inbounds @fastmath function random_three_cycle!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)    
    
    if !reverse
        # Give random cubelet_subsystem_label, and random cubelet indices (within number_of_cubelets_in_subsystem)
        cubelet_subsystem_label = cube.cubelet_subsystems_labels[rand(1:length(cube.cubelet_subsystems_labels))]

        if cubelet_subsystem_label=="sigma_"
            number_of_cubelets_in_subsystem = 8
        elseif cubelet_subsystem_label=="tau_"
            number_of_cubelets_in_subsystem = 12
        else
            number_of_cubelets_in_subsystem = 24
        end

        # Get three non-duplicated indices within the number of cubelets in the subsystem
        random_cubelet_indices = sample(1:number_of_cubelets_in_subsystem, 3, replace=false)

        three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[2], random_cubelet_indices[3])

        candidate_reversing_information = (cubelet_subsystem_label,random_cubelet_indices)

        return candidate_reversing_information

    # Reversing case
    else
        cubelet_subsystem_label, random_cubelet_indices = candidate_reversing_information

        # Now just want 3-cycle in opposite direction therefore just define the cubelets in an odd signature order e.g. [1,3,2]
        # i.e. before we did [1,2,3] ---> [3,1,2]
        # and now we are saying [(1),(3),(2)] = [3,2,1] --> [1,3,2] = [(2),(1),(3)]
        three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[3], random_cubelet_indices[2])

    end
end





@inline @inbounds @fastmath function random_orientation_rotation!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)    
    
    if !reverse

        # Give random cubelet_subsystem_label, and random cubelet indices (within number_of_cubelets_in_subsystem)

        if isodd(cube.L)
            rotatable_cubelet_subsystem_labels = ["sigma_", "tau_"]                
            cubelet_subsystem_label = rotatable_cubelet_subsystem_labels[rand(1:2)]
        else
            cubelet_subsystem_label = "sigma_"
        end


        if cubelet_subsystem_label=="sigma_"
            number_of_cubelets_in_subsystem = 8
        else 
            #(cubelet_subsystem_label=="tau_" case)
            number_of_cubelets_in_subsystem = 12
        end

        random_cubelet_indices = sample(1:number_of_cubelets_in_subsystem, 2, replace=false)

        opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[2])

        candidate_reversing_information = (cubelet_subsystem_label,random_cubelet_indices)

        return candidate_reversing_information

    # Reversing case
    else
        cubelet_subsystem_label, random_cubelet_indices = candidate_reversing_information

        # Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa
        opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[2], random_cubelet_indices[1])
    end
end





@inline function random_swap_move!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)

    # With (cube) arguments: performs random swap_move of type:
    # - P_3{X,(a,b,c)} i.e. a 3-cycle in cubelet subsystem X of cubelets in positions -->a-->b-->c-->
    # or
    # - O_{X,(a,b)} i.e. an opposite orientation rotation in cubelet subsystem X of factor +1 (anticlockwise) to cubelet in position a and -1(clockwise) to cubelet in position b
    # or
    # - Coupled P2_{X,(a,b)} i.e. coupled subsystem 2-cycles in such a way that fundamental constraints are conserved

    # Returns candidate_reversing_information = (X, (a,b) or (a,b,c))

    # With (cube, reverse=true, candidate_reversing_information) as arguments, it simply undoes the previous swap_move

    # Currently we have a perfect uniform distribution of suggesting each type of swap move, but we can bias some classes of swap moves over others using the variables below
    bias_of_coupled_subsystem_two_cycles_over_three_cycles_or_opposite_orientation_rotations = 0.0
    bias_of_three_cycles_over_opposite_orientation_rotations = 0.0


    if !reverse 
        # (i.e. if not reversing)
        # First choose random type of swap_move

        # Firstly choose whether doing a coupled subsystem 2-cycle swap move or not
        if rand() < (number_of_coupled_subsystem_two_cycle_swap_moves(cube; inequivalent=true)/total_number_of_swap_moves(cube;inequivalent=true)) + bias_of_coupled_subsystem_two_cycles_over_three_cycles_or_opposite_orientation_rotations
            # Parity sector exchange case

            candidate_reversing_information = random_parity_sector_switch_coupled_two_cycles!(cube)

        else
            # Not parity sector exchange case

            # Now choose between three cycles or opposite orientation rotations
            
            if rand() < (number_of_three_cycle_swap_moves(cube;inequivalent=true)/(number_of_three_cycle_swap_moves(cube;inequivalent=true)+number_of_opposite_orientation_rotation_swap_moves(cube))) + bias_of_three_cycles_over_opposite_orientation_rotations
                # P_{X,(a,b,c)} i.e. three cycle case 
                candidate_reversing_information = random_three_cycle!(cube)

            else 
                # O_{X,(a,b)} i.e. opposite orientation rotation case
                candidate_reversing_information = random_orientation_rotation!(cube)
                
            end

        end


    else # Reversing case

        # First deal with not parity exchange sector swap case

        if typeof(candidate_reversing_information) == Tuple{String, Vector{Int64}}

            cubelet_subsystem_label, random_cubelet_indices = candidate_reversing_information
        
            if length(random_cubelet_indices) == 3 # P_{X,(a,b,c)} case

                random_three_cycle!(cube; reverse=true, candidate_reversing_information=candidate_reversing_information)               

            else # O_{X,(a,b)} case

                random_orientation_rotation!(cube; reverse=true, candidate_reversing_information=candidate_reversing_information)

            end
        
        else  # Parity exchange sector swap case

            random_parity_sector_switch_coupled_two_cycles!(cube; reverse=true, candidate_reversing_information=candidate_reversing_information)

        end

    end

end 



# Numbers ---

function number_of_opposite_orientation_rotation_swap_moves(cube::RubiksCube)
    
    number_of_opposite_orientation_rotation_swap_moves = 0

    if cube.L == 1
        return 0
    end

    # Add for 'corners' (always 8 cubelets in subsystem, choose permutations of 2 from this 8 as clockwise/anticlockwise matters, then 2 types of rotation as 1 unit or 2 units)
    number_of_opposite_orientation_rotation_swap_moves += length(permutations(1:8, 2))*2

    if isodd(cube.L)
        # Add for 'central edges' for odd cubes (always 12 cubelets in subsystem, choose permutations of 2 from thi 12 as clockwise/anticlockiwse matters, then 1 type of rotation)
        number_of_opposite_orientation_rotation_swap_moves += length(permutations(1:12, 2))
    end

    return number_of_opposite_orientation_rotation_swap_moves

end



function number_of_three_cycle_swap_moves(cube::RubiksCube; inequivalent=false)
    
    number_of_three_cycle_swap_moves = 0

    for cubelet_subsystem_label in cube.cubelet_subsystems_labels

        if cubelet_subsystem_label=="sigma_"
            number_of_cubelets_in_subsystem = 8
        elseif cubelet_subsystem_label=="tau_"
            number_of_cubelets_in_subsystem = 12
        else
            number_of_cubelets_in_subsystem = 24
        end

        # For every subsystem add a three cycle using permutation of 3 from the number of cubelets in subystem (as order matters)
        # This will 'overcount' 'inequivalent' 3-cycles by a factor of 3 as (A,B,C) (C,A,B), (B,C,A) will have the same effect on the cube
        # However in reality every 3-cycle is overcounted by the same factor of 3 so we just use this to perform random three-cycle swap moves on cubes
        number_of_three_cycle_swap_moves += length(permutations(1:number_of_cubelets_in_subsystem, 3))

    
    end

    if !inequivalent

        return number_of_three_cycle_swap_moves

    else
        # However for calculating weighting factor between classes of swap moves (three cycles, opposite orientation rotations and coupled subsystem two cycles)
        # it is worth only accounting inequivalent 3 cycles so we do not unintentionally generate a bias towards them
        return Int(number_of_three_cycle_swap_moves/3)
    end

end



function number_of_coupled_subsystem_two_cycle_swap_moves(cube::RubiksCube; inequivalent=false)
    if cube.L == 1
        return 0
    end


    number_of_coupled_subsystem_two_cycle_swap_moves = 0

    # The number of subclasses of coupled subsystem two cycle swap moves is equivalent to the number of ways of assining +1/-1 to {sigma, {theta_k}}
    # (as by knowing the parity changes of these, we can uniquely determine the required parity changes of others)

    # Note that theta_k has k run in 1:ceil(Int, (cube.L-3)/2)
    number_of_coupled_subsystem_two_cycle_subclasses = 2^(1+ceil(Int, (cube.L-3)/2))


    # Now we generate a set of true/false values (equivalently a bitstring of 1s and 0s) for every one of these subclass cases
    # (but trim so only (1+ceil(Int, (cube.L-3)/2)) bits)
    all_coupled_subsystem_two_cycle_subclass_sigma_and_theta_k_parity_changes = [reverse(reverse(bitstring(a))[1:(1+ceil(Int, (cube.L-3)/2))]) for a in 0:number_of_coupled_subsystem_two_cycle_subclasses-1]


    # For every subclass of coupled subsystem two cycle swap move:
    for subclass in all_coupled_subsystem_two_cycle_subclass_sigma_and_theta_k_parity_changes

        # We first get the list of ALL subsystems that must be two cycled in order to respect the fundamental constraints of the cube
        two_cycle_sigma = parse(Bool,subclass[1])
        two_cycle_theta_ks = isempty(subclass[2:end]) ? [] : parse.(Bool, split(subclass[2:end],""))

        subsystems_to_two_cycle::Vector{String} = get_coupled_subsystems_to_two_cycle(cube,two_cycle_sigma, two_cycle_theta_ks)

        # For every subsystem that must be two cycled calculate and add permutations of 2 cubelets from number of cubelets in subsystem
        # TODO Note

        for subsystem_to_two_cycle in subsystems_to_two_cycle
            if subsystem_to_two_cycle=="sigma_"
                number_of_cubelets_in_subsystem = 8
            elseif subsystem_to_two_cycle=="tau_"
                number_of_cubelets_in_subsystem = 12
            else
                number_of_cubelets_in_subsystem = 24
            end

            # This will 'overcount' 'inequivalent' 2-cycles by a factor of 3 as (A,B) and (B,A) will have the same effect on the cube
            # However in reality every 2-cycle is overcounted by the same factor of 2 so we just use this to perform random two-cycle swap moves on cubes
            number_of_coupled_subsystem_two_cycle_swap_moves += length(permutations(1:number_of_cubelets_in_subsystem, 2))

        end

    end

    if !inequivalent

        return number_of_coupled_subsystem_two_cycle_swap_moves

    else
        # However for calculating weighting factor between classes of swap moves (three cycles, opposite orientation rotations and coupled subsystem two cycles)
        # it is worth only accounting inequivalent 2 cycles so we do not unintentionally generate a bias towards them
        return Int(number_of_coupled_subsystem_two_cycle_swap_moves/2)
    end

end


function total_number_of_swap_moves(cube::RubiksCube; inequivalent=false)

    total_number_of_swap_moves = 0

    total_number_of_swap_moves += number_of_opposite_orientation_rotation_swap_moves(cube)

    total_number_of_swap_moves += number_of_three_cycle_swap_moves(cube; inequivalent=inequivalent)

    total_number_of_swap_moves += number_of_coupled_subsystem_two_cycle_swap_moves(cube; inequivalent=inequivalent)

    return total_number_of_swap_moves

end


@inline function configuration_network_degree(L::Int64, including_swap_moves::Bool)

    if including_swap_moves
        return total_number_of_swap_moves(RubiksCube(L))
    else
        return total_number_of_slice_rotations(L)
    end

end 





# ----- Neighbours -----

function add_opposite_orientation_rotation_neighbour_energy_deltas!(cube::RubiksCube, opposite_orientation_rotation_neighbour_energy_deltas)
    
    current_cube_energy = energy(cube)
    neighbour_index = 1

    # Determine which rotable subsystems are present in the cube
    if isodd(cube.L)
        cubelet_subsystem_labels = ["sigma_", "tau_"]
    else
        cubelet_subsystem_labels = ["sigma_"]
    end

    # For every rotable subsystem
    for cubelet_subsystem_label in cubelet_subsystem_labels

        if cubelet_subsystem_label == "sigma_"
            number_of_cubelets_in_subsystem = 8
        else # (cubelet_subsystem_label=="tau_" case)
            number_of_cubelets_in_subsystem = 12
        end

        # For every pair of cubelets in the subsystem
        all_random_cubelet_indices_combinations = permutations(1:number_of_cubelets_in_subsystem, 2)

        for random_cubelet_indices_combination in all_random_cubelet_indices_combinations
            
            # Do one unit opposite orientation rotation to neighbouring configuration and add delta energy
            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2])
            opposite_orientation_rotation_neighbour_energy_deltas[neighbour_index] = energy(cube) - current_cube_energy
            neighbour_index += 1

            # Do the second unit of opposite orientation rotations that are availbe to corners only if evaluating for corners
            if cubelet_subsystem_label == "sigma_"
                opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2])
                opposite_orientation_rotation_neighbour_energy_deltas[neighbour_index] = energy(cube) - current_cube_energy
                neighbour_index += 1

                # Reverse second unit of rotation if applied
                # (Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa)
                opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[2], random_cubelet_indices_combination[1])
            end

            # Reverse first unit of rotation
            # (Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa)
            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[2], random_cubelet_indices_combination[1])
        end

    end
end



function add_three_cycle_swap_move_neighbour_energy_deltas!(cube::RubiksCube, three_cycle_swap_move_neighbour_energy_deltas)
    
    current_cube_energy = energy(cube)
    neighbour_index = 1




    # ---
    
    # Determine which rotable subsystems are present in the cube
    if isodd(cube.L)
        cubelet_subsystem_labels = ["sigma_", "tau_"]
    else
        cubelet_subsystem_labels = ["sigma_"]
    end

    # For every rotable subsystem
    for cubelet_subsystem_label in cubelet_subsystem_labels

        if cubelet_subsystem_label == "sigma_"
            number_of_cubelets_in_subsystem = 8
        else # (cubelet_subsystem_label=="tau_" case)
            number_of_cubelets_in_subsystem = 12
        end

        # For every pair of cubelets in the subsystem
        all_random_cubelet_indices_combinations = permutations(1:number_of_cubelets_in_subsystem, 2)

        for random_cubelet_indices_combination in all_random_cubelet_indices_combinations
            
            # Do one unit opposite orientation rotation to neighbouring configuration and add delta energy
            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2])
            opposite_orientation_rotation_neighbour_energy_deltas[neighbour_index] = energy(cube) - current_cube_energy
            neighbour_index += 1

            # Do the second unit of opposite orientation rotations that are availbe to corners only if evaluating for corners
            if cubelet_subsystem_label == "sigma_"
                opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[1], random_cubelet_indices_combination[2])
                opposite_orientation_rotation_neighbour_energy_deltas[neighbour_index] = energy(cube) - current_cube_energy
                neighbour_index += 1

                # Reverse second unit of rotation if applied
                # (Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa)
                opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[2], random_cubelet_indices_combination[1])
            end

            # Reverse first unit of rotation
            # (Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa)
            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices_combination[2], random_cubelet_indices_combination[1])
        end

    end
end



function swap_move_neighbour_energy_deltas!(cube::RubiksCube, neighbour_energy_deltas::Vector{Float64})


    # Firstly add opposite orientation rotation neighbour delta energies
    @views add_opposite_orientation_rotation_neighbour_energy_deltas!(cube, neighbour_energy_deltas[1:number_of_opposite_orientation_rotation_swap_moves(cube)])

    println(neighbour_energy_deltas[1:number_of_opposite_orientation_rotation_swap_moves(cube)])

    # Next add three cycle swap move neighbour delta energies
    add_three_cycle_swap_move_neighbour_energy_deltas!(cube, neighbour_energy_deltas[(number_of_opposite_orientation_rotation_swap_moves(cube)+1):(number_of_opposite_orientation_rotation_swap_moves(cube)+number_of_three_cycle_swap_moves(cube))])

    # # Finally add coupled subsystem two cycle swap move neighbour delta energies
    # add_coupled_subsystem_two_cycle_swap_move_neighbour_energy_deltas!(cube, neighbour_energy_deltas[(number_of_opposite_orientation_rotation_swap_moves(cube)+number_of_three_cycle_swap_moves(cube)+1):end])


    # ----

    # for f in 1:6
    #     for l in 0:floor(Int,(cube.L/2)-1)
    #         for o in 0:1
    #             # Do rotation
    #             rotate!(cube, f, l, o)

    #             # Calculate energy difference
    #             neighbour_energy_deltas[neighbour_index] = energy(cube) - current_cube_energy

    #             # Undo rotation and increase neighbour index
    #             rotate!(cube, f, l, mod(o+1,2))
    #             neighbour_index += 1
    #         end
    #     end
    # end

end







function neighbour_energy_deltas(cube::RubiksCube, including_swap_moves::Bool)

    neighbour_energy_deltas = zeros(configuration_network_degree(cube.L, including_swap_moves))

    if including_swap_moves
        swap_move_neighbour_energy_deltas!(cube, neighbour_energy_deltas)
    else
        slice_rotation_neighbour_energy_deltas!(cube, neighbour_energy_deltas)
    end

end


