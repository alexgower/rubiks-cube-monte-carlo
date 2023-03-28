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




# --- To Do: ---

# - Put diagrams from Draft 1 neat somewhere for conventions

# ----

using StatsBase

# TODO remove
using Combinatorics

include("rubiks_cube.jl")

function get_cubelet_subsystem(cube::RubiksCube,cubelet_subsystem_label::String)
    # Function to return an array of references to cube's configuration for the cubelets in subsystem_label

    # First chop cubelet_subsystem_label X_i,j into subsystem_type = X and subsystem_details= [i,j]
    subsystem_type = split(cubelet_subsystem_label,"_")[1]
    subsystem_details = split(split(cubelet_subsystem_label,"_")[2], ",")

    # Validation (only odd L cubes have central edge or +-centre cubelets)
    if (subsystem_type=="tau" && iseven(cube.L)) || (subsystem_type=="eta" && iseven(cube.L))
        throw(ArgumentError("Only odd L cubes have central edge or +-centre cubelets "))
    end

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
                            [view(cube.configuration[2],1,1), view(cube.configuration[3],cube.L,1), view(cube.configuration[5],1,1)],
                            [view(cube.configuration[2],1,cube.L), view(cube.configuration[6],1,cube.L), view(cube.configuration[3],cube.L,cube.L)],
                            [view(cube.configuration[4],cube.L,1), view(cube.configuration[5],cube.L,1), view(cube.configuration[3],1,1)],
                            [view(cube.configuration[4],cube.L,cube.L), view(cube.configuration[3],1,cube.L), view(cube.configuration[6],cube.L,cube.L)]]
        
    
    
    
    elseif subsystem_type == "theta"
        # Note each wing edge cubelet contains 2 facelets so each element in the subsystem_array will be a 2x1 array of 2 facelet views
        # The order of facelets in this 2x1 array will be in an aribtrary but well defined order (as shown in supplementary document)
        # Note we don't require a specific facelet order here as will only need to flip order later
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)


        i = parse(Int64,subsystem_details[1])

        subsystem_array = [[view(cube.configuration[3],cube.L,1+i), view(cube.configuration[2],1,1+i)],
                            [view(cube.configuration[3],cube.L,cube.L-i), view(cube.configuration[2],1,cube.L-i)],
                            [view(cube.configuration[2],cube.L,1+i), view(cube.configuration[1],1,1+i)],
                            [view(cube.configuration[2],cube.L,cube.L-i), view(cube.configuration[1],1,cube.L-i)],
                            [view(cube.configuration[1],cube.L,1+i), view(cube.configuration[4],1,1+i)],
                            [view(cube.configuration[1],cube.L,cube.L-i), view(cube.configuration[4],1,cube.L-i)],
                            [view(cube.configuration[4],cube.L,1+i), view(cube.configuration[3],1,1+i)],
                            [view(cube.configuration[4],cube.L,cube.L-i), view(cube.configuration[3],1,cube.L-i)],
                            [view(cube.configuration[5],1+i,cube.L), view(cube.configuration[1],1+i,1)],
                            [view(cube.configuration[5],cube.L-i,cube.L), view(cube.configuration[1],cube.L-i,1)],
                            [view(cube.configuration[1],1+i,cube.L), view(cube.configuration[6],1+i,1)],
                            [view(cube.configuration[1],cube.L-i,cube.L), view(cube.configuration[6],cube.L-i,1)],
                            [view(cube.configuration[5],1,1+i), view(cube.configuration[2],1+i,1)],
                            [view(cube.configuration[5],1,cube.L-i), view(cube.configuration[2],cube.L-i,1)],
                            [view(cube.configuration[6],1,1+i), view(cube.configuration[2],cube.L-i,cube.L)],
                            [view(cube.configuration[6],1,cube.L-i), view(cube.configuration[2],1+i,cube.L)],
                            [view(cube.configuration[4],1+i,cube.L), view(cube.configuration[6],cube.L,1+i)],
                            [view(cube.configuration[4],cube.L-i,cube.L), view(cube.configuration[6],cube.L,cube.L-i)],
                            [view(cube.configuration[4],1+i,1), view(cube.configuration[5],cube.L,cube.L-i)],
                            [view(cube.configuration[4],cube.L-i,1), view(cube.configuration[5],cube.L,1+i)],
                            [view(cube.configuration[6],1+i,cube.L), view(cube.configuration[3],cube.L-i,cube.L)],
                            [view(cube.configuration[6],cube.L-i,cube.L), view(cube.configuration[3],1+i,cube.L)],
                            [view(cube.configuration[5],1+i,1), view(cube.configuration[3],cube.L-i,1)],
                            [view(cube.configuration[5],cube.L-i,1), view(cube.configuration[3],1+i,1)]]
    
    elseif subsystem_type == "tau"
        # Note each central edge cubelet contains 2 facelets so each element in the subsystem_array will be a 2x1 array of 2 facelet views
        # The order of facelets in this 2x1 array will be in an aribtrary but well defined order (as shown in supplementary document)
        # Note we don't require a specific facelet order here as only flips are possible anyway
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        subsystem_array = [[view(cube.configuration[1],1,Int((cube.L+1)/2)), view(cube.configuration[2],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[1],Int((cube.L+1)/2),cube.L), view(cube.configuration[6],Int((cube.L+1)/2),1)],
                            [view(cube.configuration[1],cube.L,Int((cube.L+1)/2)), view(cube.configuration[4],1,Int((cube.L+1)/2))],
                            [view(cube.configuration[1],Int((cube.L+1)/2),1), view(cube.configuration[5],Int((cube.L+1)/2),cube.L)],
                            [view(cube.configuration[2],1,Int((cube.L+1)/2)), view(cube.configuration[3],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[3],1,Int((cube.L+1)/2)), view(cube.configuration[4],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[6],1,Int((cube.L+1)/2)), view(cube.configuration[2],Int((cube.L+1)/2),cube.L)],
                            [view(cube.configuration[6],Int((cube.L+1)/2),cube.L), view(cube.configuration[3],Int((cube.L+1)/2),cube.L)],
                            [view(cube.configuration[5],1,Int((cube.L+1)/2)), view(cube.configuration[2],Int((cube.L+1)/2),1)],
                            [view(cube.configuration[5],Int((cube.L+1)/2),1), view(cube.configuration[3],Int((cube.L+1)/2),1)],
                            [view(cube.configuration[4],Int((cube.L+1)/2),1), view(cube.configuration[5],cube.L,Int((cube.L+1)/2))],
                            [view(cube.configuration[4],Int((cube.L+1)/2),cube.L), view(cube.configuration[6],cube.L,Int((cube.L+1)/2))]]    
    
    elseif subsystem_type == "x"
        # Note each X-centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        i = parse(Int64,subsystem_details[1])

        for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,1+i)], [view(cube.configuration[face_number],1+i,cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,cube.L-i)], [view(cube.configuration[face_number],cube.L-i,1+i)]])
        end


    elseif subsystem_type == "omega"
        # Note each oblique centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)
    
        i = parse(Int64,subsystem_details[1])
        j = parse(Int64,subsystem_details[2])

        for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,1+j)], [view(cube.configuration[face_number],1+j,cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,cube.L-j)], [view(cube.configuration[face_number],cube.L-j,1+i)]])
        end



    elseif subsystem_type == "eta"
        # Note each +-centre cubelet is just a single facelet so each element in the subsystem_array will be a 1x1 array of 1 facelet views
        # The order of the cubelets in this overall subsystem_array is arbitrary but well defined (as shown in supplementary document)

        i = parse(Int64,subsystem_details[1])

        for face_number in 1:6
            subsystem_array = vcat(subsystem_array, [[view(cube.configuration[face_number],1+i,Int((cube.L+1)/2))], [view(cube.configuration[face_number],Int((cube.L+1)/2),cube.L-i)],
                                                    [view(cube.configuration[face_number],cube.L-i,Int((cube.L+1)/2))], [view(cube.configuration[face_number],Int((cube.L+1)/2),1+i)]])
        end



    end
    
    return subsystem_array


end



function three_cycle_cubelets!(cube::RubiksCube, cubelet_subsystem_label::String, cubelet_index_1::Int64, cubelet_index_2::Int64, cubelet_index_3::Int64)
    # Conducts P_{X,(a,b,c)} i.e. a 3-cycle in cubelet subsystem X of cubelets in positions -->1-->2-->3-->

    # Validation (all cubelet subsystems have size 24 except for the corners (8) and central edges (12))
    if max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 24 || min(cubelet_index_1, cubelet_index_2, cubelet_index_3) < 1 || (cubelet_subsystem_label=="sigma_" && max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 8 ) || (cubelet_subsystem_label=="tau_" && max(cubelet_index_1, cubelet_index_2, cubelet_index_3) > 12 )
        throw(ArgumentError("Cubelet indices are outside range of cubelet subsystem"))
    end

    cubelet_subsystem = get_cubelet_subsystem(cube, cubelet_subsystem_label)

    # TODO: is there a neater way to do this weird mutation using pointers in julia?
    # No matter how many facelets each cubelet in this cubelet_subsystem contains (hence the for loop):
    # Move cubelet 2's value to cubelet 3's position 
    # Move cubelet 1's value to cubelet 2's position 
    # Move cubelet 3's value to cubelet 1's position 
    for facelet_index in eachindex(cubelet_subsystem[cubelet_index_1])
        temp = copy(cubelet_subsystem[cubelet_index_3][facelet_index])

        cubelet_subsystem[cubelet_index_3][facelet_index] .= cubelet_subsystem[cubelet_index_2][facelet_index]
        cubelet_subsystem[cubelet_index_2][facelet_index] .= cubelet_subsystem[cubelet_index_1][facelet_index]
        cubelet_subsystem[cubelet_index_1][facelet_index] .= temp
    end


end



function opposite_rotate_cubelets!(cube::RubiksCube, cubelet_subsystem_label::String, cubelet_index_1, cubelet_index_2)
    # Conducts O_{X,(a,b),o} i.e. an orientation rotation in cubelet subsystem X of factor +1 (anticlockwise) to cubelet in position 1 and -1 (clockwise) to cubelet in position 2
    # e.g. [facelet_1, facelet_2, facelet_3] --> [facelet_3, facelet_1, facelet_2] is 1 (anticlockwise) rotation unit (as we used anticlockwise convention in defining facelets in cubelet array)
    # Note 

    # Note wing edge orientations are actually fully constrained by their positions so we only have corner and central edge cases that can have these orientation swap moves
    # Note corners have 3 facelets therefore can be rotated clockwise or anticlockwise, but switching cubelet_index_1 and cubelet_index_2 allows both possibilities
    # Note central edgees only have 2 facelets therefore all 'rotations' are equivalent to just flipping the facelets within the cubelet, and switching cubelet_index_1 and cubelet_index_2 makes no difference

    # Validation:
    # Only allow corners or central edges as subystems in this function, and only allow central edges for odd L cubes (as don't exist for even L cubes)
    if !(cubelet_subsystem_label=="sigma_" || cubelet_subsystem_label=="tau_") || (cubelet_subsystem_label=="tau_" && iseven(cube.L))
        throw(ArgumentError("This cubelet subsystem label ($cubelet_subsystem_label) doesn't make sense for this cube, or doesn't make sense to have an opposite rotation swap move"))
    end
    
    cubelet_subsystem = get_cubelet_subsystem(cube,cubelet_subsystem_label)

    # TODO is there a nicer way to do this without hardcoding?
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




function random_swap_move!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)

    # With (cube) arguments: performs random swap_move of type:
    # - P_{X,(a,b,c)} i.e. a 3-cycle in cubelet subsystem X of cubelets in positions -->a-->b-->c-->
    # or
    # - O_{X,(a,b)} i.e. an orientation rotation in cubelet subsystem X of factor +1 (anticlockwise) to cubelet in position a and -1(clockwise) to cubelet in position b

    # Returns candidate_reversing_information = (X, (a,b) or (a,b,c))

    # With (cube, reverse=true, candidate_reversing_information) as arguments, it simply undoes the previous swap_move

    if !reverse # (i.e. if not reversing)

        # First choose random type of swap_move

        # Make the 3-cycle swap moves equally likely to the orientation swap moves
        # (i.e uniform random choice over 3-cycle swap moves in any of the cubelet_subystems or of the 2 (corners or central edges) opposite-rotation swap moves)
        # (This isn't perfect - as not accounting for number of distinct 2/3 cubelets can pick in each cubelet_subystem for a 3-cycle/opposite-rotaiton swap move)
        # (But is much better than having a 50/50 chance between either one of the (2) types of opposite-rotation swap move, or one of the (many) types of 3-cycle swap move)
        swap_move_type = rand(1:2 + length(cube.cubelet_subsystems_labels))

        if swap_move_type > 2 # P_{X,(a,b,c)} case
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


        else # O_{X,(a,b)} case
            # Give random cubelet_subsystem_label, and random cubelet indices (within number_of_cubelets_in_subsystem)

            if isodd(cube.L)
                rotatable_cubelet_subsystem_labels = ["sigma_", "tau_"]
                cubelet_subsystem_label = rotatable_cubelet_subsystem_labels[rand(1:2)]
            else
                cubelet_subsystem_label = "sigma_"
            end


            if cubelet_subsystem_label=="sigma_"
                number_of_cubelets_in_subsystem = 8
            else #(cubelet_subsystem_label=="tau_" case)
                number_of_cubelets_in_subsystem = 12
            end

            random_cubelet_indices = sample(1:number_of_cubelets_in_subsystem, 2, replace=false)

            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[2])

        end

        return (cubelet_subsystem_label,random_cubelet_indices)

    else # Reversing case
        cubelet_subsystem_label, random_cubelet_indices = candidate_reversing_information
        
        if length(random_cubelet_indices) == 3 # P_{X,(a,b,c)} case
            # TODO remove
            # # Now just want 3-cycle in opposite direction therefore just 3-cycle twice (referencing cubelet indices in same order)
            # # i.e. before we did [1,2,3] ---> [3,1,2]
            # # and now we do [3,1,2] --> [2,3,1] --> [1,2,3] to return
            # three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[2], random_cubelet_indices[3])
            # three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[2], random_cubelet_indices[3])
            
            
            # Now just want 3-cycle in opposite direction therefore just define the cubelets in an odd signature order e.g. [1,3,2]
            # i.e. before we did [1,2,3] ---> [3,1,2]
            # and now we are saying [(1),(3),(2)] = [3,2,1] --> [1,3,2] = [(2),(1),(3)]
            three_cycle_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[1], random_cubelet_indices[3], random_cubelet_indices[2])

        else # O_{X,(a,b)} case
            # Just reverse cubelet_index order so previously anticlockwise rotated cubelet is rotated clockwise by same amount and vice versa
            opposite_rotate_cubelets!(cube, cubelet_subsystem_label, random_cubelet_indices[2], random_cubelet_indices[1])
        end

    end

end 