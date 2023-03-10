# --- To Do: ---

# - Check rotations on real Rubik's Cube

# - make energy function faster (strip half of if statements and remove half?)

# - make mode where only introduce swaps after cetain temperature



# ------ RubiksCube Struct and Configuration ------

mutable struct RubiksCube
    # An L x L x L Rubik's Cube type has a: size and configuration
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

    # Constructor 
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



function face(cube::RubiksCube, face_number::Int64)
    # Function to get an individual face_number's face configuration from a RubiksCube
    return cube.configuration[face_number]
end

function configuration_correlation_function(cube::RubiksCube, reference_cube::RubiksCube)
    # Calculates configuration correlation function with current Rubik's cube configuration compared to a cube in a reference configuration
    # The configuration correlation function is given by:
    # C_{\mathcal{C}} = \frac{1}{N_{facelets}} \sum_{facelets} \delta_{\mathcal{C}(t=0)[facelet], \mathcal{C}(t=t)[facelet]}

    # Validation - throw an error if the reference cube is not the same size
    if reference_cube.L != cube.L
        throw(ArgumentError("Reference cube has different size to current cube"))
    end

    # Calculate configuration correlation function
    number_facelets = 6*cube.L^2

    unnormalised_configuration_correlation_function = 0

    # Sum over facelets on cube
    for (face_number, face) in pairs(cube.configuration)
        for (index, current_spin) in pairs(face)
            # Implement Kronecker Delta
            if current_spin == reference_cube.configuration[face_number][index]
                unnormalised_configuration_correlation_function += 1
            end
        end
    end

    return (1/number_facelets)*unnormalised_configuration_correlation_function
end





# ------ Energy ------

# Function to get the energy of a RubiksCube 
function energy(cube::RubiksCube)

    # Start energy at 0 then count all the bonds (same spin values for nearest neighbours) present on the cube
    # with its current configuration.
    E = 0

    # Sum energy over every face individually (i.e. no couplings around corners/edges)
    for face_number in 1:6

        # Consider every facelet in the face at site (i,j) on the face (with (i,j) indexed like matrices)
        for i in 1:cube.L
            for j in 1:cube.L

                # Add up nearest neighbour (if exists) coupling energy (=-1 if same spin as current site, else 0)
                if (i-1) >= 1
                    if face(cube, face_number)[i-1, j] == face(cube,face_number)[i, j]
                        E -= 1
                    end
                end

                # Add down nearest neighbour (if exists) coupling energy (=-1 if same spin as current site, else 0)
                if (i+1) <= cube.L
                    if face(cube, face_number)[i+1, j] == face(cube, face_number)[i, j]
                        E -= 1
                    end
                end

                # Add left nearest neighbour (if exists) coupling energy (=-1 if same spin as current site, else 0)
                if (j-1) >= 1
                    if face(cube, face_number)[i, j-1] == face(cube, face_number)[i, j]
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

    # Finally return total energy value for the whole Rubik's cube 
    # This is equal to HALF of this sum which double counts bonds
    return Int(0.5*E)
end

# Function (just for convenience) to get energy of solved configuration of RubiksCube of this size left
function solved_configuration_energy(cube::RubiksCube)
    return -12*cube.L*(cube.L - 1)
end

# Function (just for convenience) to get infinite temperature energy of RubiksCube of this size
# This is equal to 1/6 * the solved configuration energy as each nearest neighbour pair has a 1/6 chance of being
#a bond (i.e. having both facelets be the same spin value = colour)
function infinite_temperature_energy(cube::RubiksCube)
    return (1/6)*solved_configuration_energy(cube)
end



# ------ Order Parameter ------

# Function to get the individual complex order parameter, m_f, for a given face, f, of the current configuration of the Rubik's cube
function face_order_parameter(cube::RubiksCube, f::Int64)
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

    # f = face number of rotation (between 1 and 6) - NOTE THIS IS A DIFFERENT CONVENTION TO THE PAPER DRAFT CURRENTLY
    # l = layer number of rotation (integer between 0 and (n+1)/2)
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
        # (This also accounts for not rotating central facelets for odd-n cubes so is consistent for both even and odd cases )
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
            # [2,cube.L-l,j]-(+)->[6,i,1+l]-(-)->[4,1+l,j]-(+)->[5,i,cube.L-l]-(-)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[2][cube.L-l,:]
            lineConfiguration2 = cube.configuration[6][:,1+l]
            lineConfiguration3 = cube.configuration[4][1+l,:]
            lineConfiguration4 = cube.configuration[5][:,cube.L-l]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[6][:, 1+l] = lineConfiguration1[:]
            cube.configuration[4][1+l, :] = reverse(lineConfiguration2[:])
            cube.configuration[5][:,cube.L-l] = lineConfiguration3[:]
            cube.configuration[2][cube.L-l,:] = reverse(lineConfiguration4[:])

        elseif f == 2

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [3,cube.L-l,j]-(-)->[6,1+l,j]-(+)->[1,1+l,j]-(+)->[5,1+l,j]-(-)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[3][cube.L-l,:]
            lineConfiguration2 = cube.configuration[6][1+l,:]
            lineConfiguration3 = cube.configuration[1][1+l,:]
            lineConfiguration4 = cube.configuration[5][1+l,:]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[6][1+l,:] = reverse(lineConfiguration1[:])
            cube.configuration[1][1+l,:] = lineConfiguration2[:]
            cube.configuration[5][1+l,:] = lineConfiguration3[:]
            cube.configuration[3][cube.L-l,:]  = reverse(lineConfiguration4[:])

        elseif f == 3

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [4,cube.L-l,j]-(+)->[6,i,1+l]-(-)->[2,1+l,j]-(-)->[5,i,1+l]-(+)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[4][cube.L-l,:]
            lineConfiguration2 = cube.configuration[6][:,1+l]
            lineConfiguration3 = cube.configuration[2][1+l,:]
            lineConfiguration4 = cube.configuration[5][:,1+l]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[6][:,1+l] = lineConfiguration1[:]
            cube.configuration[2][1+l,:] = reverse(lineConfiguration2[:])
            cube.configuration[5][:,1+l] = reverse(lineConfiguration3[:])
            cube.configuration[4][cube.L-l,:] = lineConfiguration4[:]

        elseif f == 4

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [1,cube.L-l,j]-(+)->[6,cube.L-l,j]-(-)->[3,1+l,j]-(-)->[5,cube.L-l,j]-(+)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[1][cube.L-l,:]
            lineConfiguration2 = cube.configuration[6][cube.L-l,:]
            lineConfiguration3 = cube.configuration[3][1+l,:]
            lineConfiguration4 = cube.configuration[5][cube.L-l,:]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[6][cube.L-l,:] = lineConfiguration1[:]
            cube.configuration[3][1+l,:] = reverse(lineConfiguration2[:])
            cube.configuration[5][cube.L-l,:] = reverse(lineConfiguration3[:])
            cube.configuration[1][cube.L-l,:]  = lineConfiguration4[:]

        elseif f == 5

            # 4 cycle of [f, i, j] (where (-) means reverse order):
            # [3,i,1+l]-(+)->[2,i,1+l]-(+)->[1,i,1+l]-(+)->[4,i,1+l]-(+)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[3][:,1+l]
            lineConfiguration2 = cube.configuration[2][:,1+l]
            lineConfiguration3 = cube.configuration[1][:,1+l]
            lineConfiguration4 = cube.configuration[4][:,1+l]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[2][:,1+l] = lineConfiguration1[:]
            cube.configuration[1][:,1+l] = lineConfiguration2[:]
            cube.configuration[4][:,1+l] = lineConfiguration3[:]
            cube.configuration[3][:,1+l] = lineConfiguration4[:]

        elseif f == 6

            # 4 cycle of [f, i, j]:
            # [1,i,cube.L-l]-(+)->[2,i,cube.L-l]-(+)->[3,i,cube.L-l]-(+)->[4,i,cube.L-l]-(+)->

            # Gather line configurations that are going to be moved
            lineConfiguration1 = cube.configuration[1][:,cube.L-l]
            lineConfiguration2 = cube.configuration[2][:,cube.L-l]
            lineConfiguration3 = cube.configuration[3][:,cube.L-l]
            lineConfiguration4 = cube.configuration[4][:,cube.L-l]

            # Assign line configurations to correct new places
            # Ensure to reverse order if need be
            cube.configuration[2][:,cube.L-l] = lineConfiguration1[:]
            cube.configuration[3][:,cube.L-l] = lineConfiguration2[:]
            cube.configuration[4][:,cube.L-l] = lineConfiguration3[:]
            cube.configuration[1][:,cube.L-l] = lineConfiguration4[:]

        end

        # Now if we are rotating a face (i.e. l=0) then we have to rotate the facelets on the actual face too
        # np.rot90() is an anti-clockwise rotation by default so we have to do it 3 times to get an clockwise rotation
        if l == 0
            cube.configuration[f, :, :] = rotr90(cube.configuration[f, :, :])
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
        l = rand(0:floor(Int,(cube.L/2) -1))

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