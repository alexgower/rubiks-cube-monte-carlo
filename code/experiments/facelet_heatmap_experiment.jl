using Plots
using Plots.PlotMeasures

include("../core/monte_carlo.jl")
include("../core/rubiks_cube.jl")
include("../core/swap_moves.jl")

function plot_facelet_heatmaps(facelet_heatmap; accepted_candidates::Int=0, T::Float64=NaN)
    nrows = 2
    ncols = 3
    plot_list = []
    max_value = maximum([maximum(facelet_heatmap[f]) for f in 1:6])  # For consistent color scaling
    
    # Slightly reduce the size of each subplot
    subplot_size = (300, 300)
    
    for face_number in 1:6
        data = facelet_heatmap[face_number]
        # Reverse the data vertically to match visual orientation
        data_reversed = reverse(data, dims=1)
        
        p = heatmap(
            data_reversed,
            title="Face $face_number",
            colorbar_title="",
            colorbar=:right,
            clims=(0, max_value),
            aspect_ratio=:equal,
            framestyle=:box,
            xticks=false,
            yticks=false,
            titlefontsize=8,
            tickfontsize=6,
            guidefontsize=8,
            size=subplot_size
        )
        push!(plot_list, p)
    end
    
    # Combine all plots into a grid layout with reduced overall size
    combined_plot = plot(plot_list...,
                         layout=(nrows, ncols),
                         size=(700, 500),  # Reduced overall size
                         link=:both,  # Link axes across subplots
                         suptitle="Heatmaps at T=$(round(T, digits=2)), Accepted Candidates: $accepted_candidates",
                         titlefontsize=10)
    
    display(combined_plot)
end

function facelet_heatmap_metropolis_swap_algorithm!(cube::RubiksCube, beta::Float64;
    swap_move_probability::Float64=0.0,
    maximum_iterations::Int64=1000,
    configuration_correlation_convergence_criteria::Float64=exp(-1),
    verbose::Bool=false,
    track_facelet_heatmap::Bool=true)

    if !track_facelet_heatmap
        accepted_candidates, _ = run_metropolis_swap_algorithm!(cube, beta;
            swap_move_probability=swap_move_probability,
            maximum_iterations=maximum_iterations,
            configuration_correlation_convergence_criteria=configuration_correlation_convergence_criteria,
            verbose=verbose)
        return (accepted_candidates, nothing)
    end

    # Initialize variables
    current_iteration = 0
    accepted_candidates = 0
    facelet_heatmap = [zeros(cube.L, cube.L) for _ in 1:6]

    # Initial configuration
    initial_cube = RubiksCube(cube.L)
    initial_cube.configuration = deepcopy(cube.configuration)
    current_configuration_correlation_function_value = configuration_correlation_function(cube, initial_cube)

    while (current_iteration <= maximum_iterations) &&
          (current_configuration_correlation_function_value > configuration_correlation_convergence_criteria)
        # Determine candidate generating function
        candidate_generating_function! = rand() < swap_move_probability ? random_swap_move! : random_rotate!

        # Perform Monte Carlo timestep
        accepted_candidates_increase, _ = monte_carlo_timestep!(cube, candidate_generating_function!, beta; verbose=verbose)

        # Update facelet heatmap
        for f in 1:6
            for i in 1:cube.L
                for j in 1:cube.L
                    if cube.configuration[f][i, j] != initial_cube.configuration[f][i, j]
                        facelet_heatmap[f][i, j] += 1
                    end
                end
            end
        end

        # Update iteration and accepted candidates
        current_iteration += 1
        accepted_candidates += accepted_candidates_increase
        current_configuration_correlation_function_value = configuration_correlation_function(cube, initial_cube)
    end

    return accepted_candidates, facelet_heatmap
end

function heatmap_anneal!(cube::RubiksCube, temperature_vector::Vector{Float64},
                         heatmap_temperatures::Vector{Float64}; anneal_swap_move_probability::Float64=1.0)
    # Initialize a dictionary to store the heatmaps
    heatmaps = Dict{Float64, Any}()

    # Mixing Stage: Randomize the cube at infinite temperature
    facelet_heatmap_metropolis_swap_algorithm!(cube, 0.0;
        swap_move_probability=0.0,
        maximum_iterations=1000,
        configuration_correlation_convergence_criteria=exp(-10),
        track_facelet_heatmap=false)

    # Annealing Stage: Loop over the temperatures
    for T in temperature_vector
        beta = 1.0 / T
        track_facelet_heatmap = T in heatmap_temperatures
        swap_move_probability = T in heatmap_temperatures ? 0.0 : anneal_swap_move_probability

        println("Running facelet_heatmap_metropolis_swap_algorithm! at temperature T = $T")

        # Run the algorithm
        accepted_candidates, facelet_heatmap = facelet_heatmap_metropolis_swap_algorithm!(cube, beta;
            swap_move_probability=swap_move_probability,
            maximum_iterations=1000,
            configuration_correlation_convergence_criteria=0.0,
            verbose=false,
            track_facelet_heatmap=track_facelet_heatmap)

        # If we tracked the heatmap, store it along with accepted candidates
        if track_facelet_heatmap
            heatmaps[T] = (accepted_candidates, facelet_heatmap)
        end
    end

    # Normalize heatmaps
    for T in heatmap_temperatures
        accepted_candidates, facelet_heatmap = heatmaps[T]
        # Normalize and round the heatmaps
        facelet_heatmap = [round.(heatmap ./ sum(heatmap); digits=2) for heatmap in facelet_heatmap]
        # Update the heatmaps dictionary
        heatmaps[T] = (accepted_candidates, facelet_heatmap)
    end

    # After collecting the heatmaps
    for T in sort(collect(keys(heatmaps)), rev=true)
        println("Heatmap at Temperature T = $T")
        accepted_candidates, facelet_heatmap = heatmaps[T]
        println("Accepted Candidates: $accepted_candidates")
        for face_number in 1:6
            println("Face $face_number Heatmap:")
            for row in eachrow(facelet_heatmap[face_number])
                println(join(row, " "))
            end
            println()  # Add an empty line between faces
        end

        plot_facelet_heatmaps(facelet_heatmap; accepted_candidates=accepted_candidates, T=T)
        savefig("heatmap_T_$(round(T, digits=2)).png")
    end

    return heatmaps
end

function run_heatmap_anneal(L::Int; T1::Float64, T0::Float64, N_T::Int64,
                            heatmap_temperatures::Vector{Float64}, anneal_swap_move_probability::Float64)
    # Initialize the Rubik's cube
    cube = RubiksCube(L)

    # Create the logarithmic decaying temperature vector
    temperature_vector = [T1 * (T0 / T1)^(m / N_T) for m in 0:N_T]
    temperature_vector = sort(vcat(temperature_vector, heatmap_temperatures), rev=true)

    # Run the heatmap annealing process
    heatmaps = heatmap_anneal!(cube, temperature_vector, heatmap_temperatures;
        anneal_swap_move_probability=anneal_swap_move_probability)

    return heatmaps
end
