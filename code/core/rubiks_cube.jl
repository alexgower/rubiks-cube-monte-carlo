using Random
using LinearAlgebra

# ------ RubiksCube Struct and Configuration ------

mutable struct RubiksCube
    # An L x L x L Rubik's Cube type has a size and configuration
    # (these are initialised together so won't be inconsistent unless reassign L maliciously)

    # Each facelet can have one of 6 spins (colours):
    #     1 = Green
    #     2 = White
    #     3 = Blue
    #     4 = Yellow
    #     5 = Orange
    #     6 = Red


    # Configuraiton of Rubik's Cube (6D vector of face configuration matrices)
    # For odd L cubes - the face configurations are indexed by the above colour-number mapping for their central facelets
    # For even L cubes - the face configurations are arbitrary (there is no fixed central facelet)
    configuration::Vector{Matrix{Int64}}

    # Size (facelets per side) of Rubik's Cube
    L::Int64 

    # Labels of independent cubelet subsystems for cube of size L (e.g. "sigma", "tau", "theta_1", "theta_2")
    cubelet_subsystems_labels::Vector{String}

    # Constructor:
    # Create Rubik's cube in solved configuration by default, and ensure fields are self consistent 
    RubiksCube(L) = new(solved_configuration(L), L, get_cubelet_subsystems_labels(L))
end


function solved_configuration(L::Int64)
    # Function to define solved_configuration for a cube of size L (i.e. face 1 has all entries of 1 etc)
    return [fill(face_number,(L,L)) for face_number in 1:6]
end


function get_cubelet_subsystems_labels(L::Int64)
    # Function to return an array of labels of independent cubelet subsystems for a Rubik's Cube of size L

    # Trivially if L==1 then there are no cubelet subsystems (just 'fixed centres' which don't count)
    if L==1 
        return []
    end

    # Add cubelet subsystems that exist for both even and odd cube.L cubes
    # Namely:
    # - Corner subsystem, sigma
    # - Wing edge subsystems, theta_k for k in 1:ceil((L-3)/2)
    # - X centre subsystems, x_i for i in 1:ceil((L-3)/2) -- Note this is labelled as omega_ii in our paper's CONVENTION
    # - Oblique centre subsystems, omega_i,j fo4 i,j in 1:ceil((L-3)/2) and i!=j
    cubelet_subsystems = vcat(["sigma_"],["theta_$k" for k in 1:ceil(Int, (L-3)/2)], ["x_$i" for i in 1:ceil(Int, (L-3)/2)], ["omega_$i,$j" for i in 1:ceil(Int, (L-3)/2) for j in 1:ceil(Int, (L-3)/2) if i!=j])

    # Add cubelet subsystems that only exist for odd cube.L case
    # Namely:
    # - Central edge subsystem, tau
    # - +-centre subsystems eta_k for k in 1:ceil((L-3)/2)
    if isodd(L)
        cubelet_subsystems = vcat(cubelet_subsystems, ["tau_"], ["eta_$k" for k in 1:ceil(Int, (L-3)/2)])
    end

    return cubelet_subsystems

end


@inline function face(cube::RubiksCube, face_number::Int64)
    # Function to get an individual face_number's face configuration from a RubiksCube
    return cube.configuration[face_number]
end

# Function with cube arguments, just pass configurations to separate function
function configuration_correlation_function(cube::RubiksCube, reference_cube::RubiksCube)
    return configuration_correlation_function(cube.configuration, reference_cube.configuration)
end

# Separate function with just the configurations as arguments
function configuration_correlation_function(cube_configuration::Vector{Matrix{Int64}}, reference_cube_configuration::Vector{Matrix{Int64}})
    # Calculates configuration correlation function with current Rubik's cube configuration compared to a cube in a reference configuration
    # The configuration correlation function is given by:
    # C_{\mathcal{C}} = \frac{1}{N_{facelets}} \sum_{facelets} \delta_{\mathcal{C}(t=0)[facelet], \mathcal{C}(t=t)[facelet]}

    # Validation ---

    # Throw an error if the reference cube is not the same size
    if size(cube_configuration) != size(reference_cube_configuration)
        throw(ArgumentError("Reference cube has different size to current cube"))
    end

    L = size(cube_configuration[1])[1]

    # Calculation ---

    # Calculate configuration correlation function
    number_facelets = 6*L^2

    unnormalised_configuration_correlation_function = 0

    # Sum over facelets on cube
    for (face_number, face) in pairs(cube_configuration)
        for (index, current_spin) in pairs(face)
            # Implement Kronecker Delta
            if current_spin == reference_cube_configuration[face_number][index]
                unnormalised_configuration_correlation_function += 1
            end
        end
    end

    # If cube is odd then we subtract 6 from the number of facelets as the central facelets are fixed
    if isodd(L)
        unnormalised_configuration_correlation_function -= 6
        number_facelets -= 6
    end

    return (1/number_facelets)*unnormalised_configuration_correlation_function
end




# ------ Energy ------

# Function to get the energy of a RubiksCube 
function energy(cube::RubiksCube)::Float64
    # Start energy at 0 then count all the bonds (same spin values for nearest neighbours) present on the cube
    # with its current configuration.
    E = 0.0

    # Sum energy over every face individually (i.e. no couplings around corners/edges)
    for face_number in 1:6

        # Consider every facelet in the face at site (i,j) on the face (with (i,j) indexed like matrices)
        for i in 1:cube.L
            for j in 1:cube.L

                # Note can remove up and left nearest neighbour bond checks as will have already been included once 
                # from another facelet (hence we don't need the factor of 1/2 in energy at the end)

                # Add down nearest neighbour (if exists) coupling energy (=-1 if same spin as current site, else 0)
                if (i+1) <= cube.L
                    if face(cube, face_number)[i+1, j] == face(cube, face_number)[i, j]
                        E -= 1
                    end
                end

                # Add right nearest neighbour (if exists) coupling energy (=-1 if same spin as current site, else 0)
                if (j+1) <= cube.L
                    if face(cube, face_number)[i, j+1] == face(cube, face_number)[i, j]
                        E -= 1
                    end
                end

            end
        end
    end

    return E
end

# Function (just for convenience) to get energy of solved configuration of RubiksCube of this size left
@inline function solved_configuration_energy(cube::RubiksCube)::Float64
    return -12*cube.L*(cube.L - 1)
end

@inline function solved_configuration_energy(L::Int64)::Float64
    return solved_configuration_energy(RubiksCube(L))
end

# Function (just for convenience) to get infinite temperature energy of RubiksCube of this size
# This is equal to 1/6 * the solved configuration energy as each nearest neighbour pair has a 1/6 chance of being
# a bond (i.e. having both facelets be the same spin value = colour)
@inline function infinite_temperature_energy(cube::RubiksCube)
    return (1/6)*solved_configuration_energy(cube)
end





# ------ Order Parameter ------

# Function to get the individual complex order parameter, m_f, for a given face, f, of the current configuration of the Rubik's cube
@inline function face_order_parameter(cube::RubiksCube, f::Int64)
    # Start face order parameter as 0 then sum over facelets on face
    m_f = 0

    # Sum over facelets on face
    for spin in face(cube, f)
        m_f += (1/cube.L)^2 * exp(((2*pi*im)/6)*spin)
    end

    return m_f
end

# Function to get the overall order parameter, M^2, of the whole RubiksCube
function order_parameter(cube::RubiksCube)
    # Start with order parameter as 0 then sum over all face order parameters
    M_squared = 0

    # Sum over all face order parameters
    for f in 1:6
        M_squared += (1/6)*abs(face_order_parameter(cube, f))^2
    end

    return M_squared
end





# ------ Rotations ------

function rotate!(cube::RubiksCube, f::Int64, l::Int64, o::Int64)
    # Performs R_{f,l,o} rotation, where the indices f and l correspond to the face and layer of the rotation (where l=0 is the face itself)
    # Note this is a mutating function i.e. it modifies the cube object in place and does not return anything

    # f = face number of rotation (between 1 and 6)
    # l = layer number of rotation (integer between 0 and (n-1)/2)
    # o = orientation number of rotation (0 = clockwise, 1 = anticlockwise)

    # Firstly notice that an anticlockwise (o=1) rotation is always equal to 3 successive clockwise (o=0) rotations
    # therefore we can implement this recursively
    # Note that when we are 'in' a recursion instance, o will equal 0 so this bit will be skipped (so no infinite recursion cycles)
    if o == 1
        # For anticlockwise (o=1) rotations

        rotate!(cube, f, l, 0)
        rotate!(cube, f, l, 0)
        rotate!(cube, f, l, 0)

    else

        # From now on only implementing clockwise rotations

        # If L=1 then we have a single odd-n cubelet which is not allowed to rotate so do nothing
        if cube.L == 1
            return
        end

        # Rotations only make sense if l+1 (accounting for 0 index for l) is less than (n+1)/2 (i.e. the middle layer)
        # (This also accounts for not rotating central facelets for odd-n cubes so is consistent for both even and odd cases)
        # Also l must be positive and f must be =< 6 and f must be >=1 and o must be 0 or 1
        if !(l+1 < (cube.L + 1)/2) || (l < 0) || (f > 6) || (f < 1) || !(o == 0 || o == 1)
            throw(ArgumentError("These rotation indices f,l,o do not make sense"))
        end

        # Now we hardcode how the rows/columns of each face configuration array are rotated into rows/columns of
        # other face configuration arrays for each of the 6 types/faces of clockwise rotation
        # (layers just change relative row/colum indices by +-l)

        # Firstly we deal with the facelet changes on the 'side' of a layer (this is all there is to a rotation
        # unless we are rotating a face itself i.e. unless l=0)

        if f == 1

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [2,cube.L-l,:] -(+)-> [6,:,1+l] -(-)-> [4,1+l,:] -(+)-> [5,:,cube.L-l] -(-)-> 

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[5][:,cube.L-l]

            cube.configuration[5][:,cube.L-l] = copy(cube.configuration[4][1+l,:])
            cube.configuration[4][1+l, :] = reverse(copy(cube.configuration[6][:,1+l]))
            cube.configuration[6][:, 1+l] = copy(cube.configuration[2][cube.L-l,:])
            cube.configuration[2][cube.L-l,:] = reverse(copy(temp))

        elseif f == 2

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [3,cube.L-l,:] -(-)-> [6,1+l,:] -(+)-> [1,1+l,:] -(+)-> [5,1+l,:] -(-)->

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[5][1+l,:]

            cube.configuration[5][1+l,:] = copy(cube.configuration[1][1+l,:])
            cube.configuration[1][1+l,:] = copy(cube.configuration[6][1+l,:])
            cube.configuration[6][1+l,:] = reverse(copy(cube.configuration[3][cube.L-l,:]))
            cube.configuration[3][cube.L-l,:] = reverse(copy(temp))

        elseif f == 3

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [4,cube.L-l,:] -(-)-> [6,:,cube.L-l] -(+)-> [2,1+l,:] -(-)-> [5,:,1+l] -(+)->

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[5][:,1+l]

            cube.configuration[5][:,1+l] = reverse(copy(cube.configuration[2][1+l,:]))
            cube.configuration[2][1+l,:] = copy(cube.configuration[6][:,cube.L-l])
            cube.configuration[6][:,cube.L-l] = reverse(copy(cube.configuration[4][cube.L-l,:]))
            cube.configuration[4][cube.L-l,:] = copy(temp)

        elseif f == 4

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [1,cube.L-l,:] -(+)-> [6,cube.L-l,:] -(-)-> [3,1+l,:] -(-)-> [5,cube.L-l,:] -(+)->

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[5][cube.L-l,:]

            cube.configuration[5][cube.L-l,:] = reverse(copy(cube.configuration[3][1+l,:]))
            cube.configuration[3][1+l,:] = reverse(copy(cube.configuration[6][cube.L-l,:]))
            cube.configuration[6][cube.L-l,:] = copy(cube.configuration[1][cube.L-l,:])
            cube.configuration[1][cube.L-l,:] = copy(temp)

        elseif f == 5

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [3,:,1+l] -(+)-> [2,:,1+l] -(+)-> [1,:,1+l] -(+)-> [4,:,1+l]-(+)->

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[4][:,1+l]

            cube.configuration[4][:,1+l] = copy(cube.configuration[1][:,1+l])
            cube.configuration[1][:,1+l] = copy(cube.configuration[2][:,1+l])
            cube.configuration[2][:,1+l] = copy(cube.configuration[3][:,1+l])
            cube.configuration[3][:,1+l] = copy(temp)

        elseif f == 6

            # 4 cycle of [f, i, j]:
            # [1,:,cube.L-l] -(+)-> [2,:,cube.L-l] -(+)-> [3,:,cube.L-l] -(+)-> [4,:,cube.L-l]-(+)->

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            temp = cube.configuration[4][:,cube.L-l]

            cube.configuration[4][:,cube.L-l] = copy(cube.configuration[3][:,cube.L-l])
            cube.configuration[3][:,cube.L-l] = copy(cube.configuration[2][:,cube.L-l])
            cube.configuration[2][:,cube.L-l] = copy(cube.configuration[1][:,cube.L-l])
            cube.configuration[1][:,cube.L-l] = copy(temp)

        end

        # Now if we are rotating a face (i.e. l=0) then we have to rotate the facelets on the actual face too
        # rot90() is a clockwise rotation by default
        if l == 0
            cube.configuration[f] .= rotr90(cube.configuration[f])
        end
    end
end

function random_rotate!(cube::RubiksCube; reverse::Bool=false, candidate_reversing_information=nothing)
    # With (cube) arguments: performs random rotation type R_{f,l,o} to the Rubik's cube, returns (f,l,o) of that rotation
    # (i.e the candidate_reversing_information)

    # With (cube, reverse=true, (f,l,o)) as arguments, it simply undoes the previous rotation, i.e. calls rotate!(cube, f, l, mod(o+1,2))

    # f = face number of rotation (between 1 and 6)
    # l = layer number of rotation (integer between 0 and (n+1)/2)
    # o = orientation number of rotation (0 = clockwise, 1 = anticlockwise)

    if !reverse
        # First choose random face to rotate
        f = rand(1:6)

        # Next choose random layer to rotate
        # We have that 0 < l < floor[n/2 - 1]  i.e. odd-n cubes cannot rotate central facelet, and even-n cubes have
        # layers up to 'half way' (after which layers are associated with opposite face)
        l = rand(0:floor(Int,(cube.L/2)-1))

        # Next choose a random orientation to rotate (0 = clockwise, 1 = anticlockwise)
        o = rand(0:1)

        # Do rotation
        rotate!(cube, f, l, o)

        return (f,l,o)
    else
        f, l, o = candidate_reversing_information
        # Note mod(o+1,2) reverses the rotation orientation direction
        rotate!(cube, f, l, mod(o+1,2))
    end
end





# ----- Numbers -----
function total_number_of_slice_rotations(L::Int64)
    # Returns the total number of swap moves in a cube of size L
    # Note that this is the same as the number of edges in the configuration network
    # (i.e. the number of edges in the configuration graph)

    if L == 1
        return 0
    end

    if isodd(L)
        return 6*(L-1)
    else
        return 6*L
    end
end

@inline function get_number_of_cubelets_in_subsystem(subsystem_name::String)

    if subsystem_name=="sigma_"
        return 8
    elseif subsystem_name=="tau_"
        return 12
    else
        return 24
    end

end





# ----- Neighbours -----

# TODO delete once incorporate normal neighbour sampling into newer function
# function all_slice_rotation_neighbour_energies!(cube::RubiksCube, neighbour_energies; recursive_additional_neighbour_steps::Int64=0, neighbour_index::Int64=1, excluded_slice_rotation=(0,0,0))

#     # Do additional neighbour steps if required
#     if recursive_additional_neighbour_steps > 0

#         # For all slice rotation neighbours
#         for f in 1:6
#             for l in 0:floor(Int,(cube.L/2)-1)
#                 for o in 0:1

#                     # Do not include excluded slice rotation (which will correspond to immediately undo-ing the previous slice rotation when we are in a recursive call)
#                     if (f,l,o) != excluded_slice_rotation


#                             # Do rotation
#                             rotate!(cube, f, l, o)

#                             # Recurse
#                             neighbour_index = all_slice_rotation_neighbour_energies!(deepcopy(cube), neighbour_energies; recursive_additional_neighbour_steps=recursive_additional_neighbour_steps-1, neighbour_index=neighbour_index, excluded_slice_rotation=(f,l,mod(o+1,2)))

#                             # Undo rotation
#                             rotate!(cube, f, l, mod(o+1,2))

#                     end
#                 end
#             end
#         end

#         return neighbour_index

#     # Else store energy and increment neighbour index
#     elseif recursive_additional_neighbour_steps==0

#         neighbour_energies[neighbour_index] = energy(cube)
#         return neighbour_index + 1
    
#     end

# end

function all_slice_rotation_energy_connections!(cube::RubiksCube, energy_connections; recursive_neighbour_order_to_measure_to::Int64=1, connection_index::Int64=1, excluded_slice_rotation=(0,0,0), verbose::Bool=false)

    # Do additional neighbour steps if required
    if recursive_neighbour_order_to_measure_to > 0

        # For all slice rotation neighbours
        for f in 1:6
            for l in 0:floor(Int,(cube.L/2)-1)
                for o in 0:1

                    # Do not include excluded slice rotation (which will correspond to immediately undo-ing the previous slice rotation when we are in a recursive call)
                    if (f,l,o) != excluded_slice_rotation


                            # Do rotation
                            rotate!(cube, f, l, o)

                            # Recurse
                            connection_index = all_slice_rotation_energy_connections!(cube, energy_connections; recursive_neighbour_order_to_measure_to=recursive_neighbour_order_to_measure_to-1, connection_index=connection_index, excluded_slice_rotation=(f,l,mod(o+1,2)))

                            # Undo rotation
                            rotate!(cube, f, l, mod(o+1,2))

                    end
                end
            end
        end

        return connection_index

    # Else store energy and increment neighbour index
    elseif recursive_neighbour_order_to_measure_to==0

        # Get current energy of cube
        E_current = energy(cube)
        # Do excluded slice rotation to get previous energy
        rotate!(cube, excluded_slice_rotation[1], excluded_slice_rotation[2], excluded_slice_rotation[3])
        E_previous = energy(cube)
        # Rotate cube back
        rotate!(cube, excluded_slice_rotation[1], excluded_slice_rotation[2], mod(excluded_slice_rotation[3]+1,2))

        # Store energy connection
        energy_connections[connection_index] = (E_previous, E_current)

        # Print energy decreasing connections if verbose
        if verbose && E_previous > E_current
            println("E_previous > E_current = ", E_previous, " > ", E_current)

            println("Current Configuration")
            println(cube.configuration)

            println("(f,l,o) = ", excluded_slice_rotation)
            println("Previous Configuration")
            rotate!(cube, excluded_slice_rotation[1], excluded_slice_rotation[2], excluded_slice_rotation[3])
            println(cube.configuration)
            
            # Now rotate back
            rotate!(cube, excluded_slice_rotation[1], excluded_slice_rotation[2], mod(excluded_slice_rotation[3]+1,2))
        end

        # Move pointer to next location in energy_connections array
        return connection_index + 1
    
    end

end


function cumulative_all_slice_rotation_energy_connections!(cube::RubiksCube, energy_connections; neighbour_order_to_measure_to::Int64=1, cumulative_energy_connection::Union{Vector{Float64},Nothing}=nothing, order_index::Int64=0, connection_index::Int64=1, excluded_slice_rotation=(0,0,0), verbose::Bool=false
)
    # If cumulative_energy_connection=nothing then set it as empty vector of length neighbour_order_to_measure_to+1
    if isnothing(cumulative_energy_connection)
        cumulative_energy_connection = zeros(neighbour_order_to_measure_to+1)
    end

    # Add this energy to the cumulative energy connection
    cumulative_energy_connection[order_index+1] = energy(cube)

    # Do additional neighbour steps if required (and keep building up the cumulate energy connection)
    if order_index < neighbour_order_to_measure_to
        # For all slice rotation neighbours
        for f in 1:6
            for l in 0:floor(Int,(cube.L/2)-1)
                for o in 0:1
                    # Do not include excluded slice rotation (which will correspond to immediately undo-ing the previous slice rotation when we are in a recursive call)
                    if (f,l,o) != excluded_slice_rotation
                        # Do rotation
                        rotate!(cube, f, l, o)

                        # Recurse
                        connection_index = cumulative_all_slice_rotation_energy_connections!(
                            cube, 
                            energy_connections; 
                            neighbour_order_to_measure_to=neighbour_order_to_measure_to, 
                            cumulative_energy_connection=cumulative_energy_connection, 
                            order_index=order_index+1, 
                            connection_index=connection_index, 
                            excluded_slice_rotation=(f,l,mod(o+1,2))
                        )

                        # Undo rotation
                        rotate!(cube, f, l, mod(o+1,2))
                    end
                end
            end
        end
    # If we've reached the order to measure to then store the complete cumulative energy connection
    else
        if connection_index <= size(energy_connections, 1)
            energy_connections[connection_index,:] .= cumulative_energy_connection
            connection_index += 1
        end
    end

    return connection_index
end