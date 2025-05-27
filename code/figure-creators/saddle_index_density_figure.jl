using Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors
using Statistics

using LsqFit

include("../core/rubiks_cube.jl")

function saddle_index_density_figure(neighbour_order_to_measure_to::Int64=1; fitting::Bool=false)

    simulation_name = "L=11_combined"

    ### -- SET UP DEFAULT PARAMETERS --
    header_line = readlines(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    energy_saddle_index_densities_data_matrix_slice = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"), ',', Float64, '\n', skipstart=3)

    # Remove any rows with energy values larger than e10
    energy_saddle_index_densities_data_matrix_slice = energy_saddle_index_densities_data_matrix_slice[energy_saddle_index_densities_data_matrix_slice[:,1] .< 1e10, :]
    # or with energy values >=0.0
    energy_saddle_index_densities_data_matrix_slice = energy_saddle_index_densities_data_matrix_slice[energy_saddle_index_densities_data_matrix_slice[:,1] .< abs(solved_configuration_energy(L)), :]
    # or with energy values <= e-100
    energy_saddle_index_densities_data_matrix_slice = energy_saddle_index_densities_data_matrix_slice[energy_saddle_index_densities_data_matrix_slice[:,1] .> -1e-10, :]
    # or with density values <= e-100
    energy_saddle_index_densities_data_matrix_slice = energy_saddle_index_densities_data_matrix_slice[energy_saddle_index_densities_data_matrix_slice[:,2] .> -1e-10, :]





    energy_saddle_index_densities_data_matrix_swap = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/",simulation_name*"_energy_saddle_index_densities_swap"), ',', Float64, '\n', skipstart=3)

    # Remove any rows with energy values larger than e10
    energy_saddle_index_densities_data_matrix_swap = energy_saddle_index_densities_data_matrix_swap[energy_saddle_index_densities_data_matrix_swap[:,1] .< 1e10, :]
    # or with energy values >=0.0
    energy_saddle_index_densities_data_matrix_swap = energy_saddle_index_densities_data_matrix_swap[energy_saddle_index_densities_data_matrix_swap[:,1] .< abs(solved_configuration_energy(L)), :]
    # or with energy values <= e-100
    energy_saddle_index_densities_data_matrix_swap = energy_saddle_index_densities_data_matrix_swap[energy_saddle_index_densities_data_matrix_swap[:,1] .> -1e-10, :]

    # E_star = -0.39015151515151514
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    # E_on = -0.23134348484848488
    E_on = -0.24010378787878783


    ### --- COLOURS ---
    Plots.default(dpi = 600)


    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    ### --- CALCULATIONS ---

    ## -- SLICE CALCULATION -- 

    # Get the unique E0 values
    E0_bin_values_slice = unique(energy_saddle_index_densities_data_matrix_slice[:,1])

    # Now do the same calculation but if we exclude K=1 saddles
    # i.e. Remove all values whose * 6(L-1) rounds to 1.0 (i.e. K=1 saddle index)
    Z = 6*(L-1)
    normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
    neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)

    # For each E0 value, find all saddle index densities and average them
    average_saddle_index_densities_slice_including_K_1_slice = zeros(Float64, length(E0_bin_values_slice))
    average_saddle_index_densities_slice_excluding_K_1_slice = zeros(Float64, length(E0_bin_values_slice))
    sem_saddle_index_densities_slice_including_K_1_slice = zeros(Float64, length(E0_bin_values_slice))
    for (i, E0) in pairs(E0_bin_values_slice)
        # Get the saddle index densities for the current E0 slice
        saddle_index_densities_including_K_1_slice = [energy_saddle_index_densities_data_matrix_slice[j,2] for j in 1:size(energy_saddle_index_densities_data_matrix_slice,1) if energy_saddle_index_densities_data_matrix_slice[j,1] == E0] 
        
        if neighbour_order_to_measure_to > 1 # Account for immediate reversal exclusion for neighbour_order_to_measure_to > 1
            saddle_index_densities_including_K_1_slice = saddle_index_densities_including_K_1_slice .+ (1/Z)
        end

        average_saddle_index_densities_slice_including_K_1_slice[i] = mean(saddle_index_densities_including_K_1_slice)
        sem_saddle_index_densities_slice_including_K_1_slice[i] = std(saddle_index_densities_including_K_1_slice) / sqrt(length(saddle_index_densities_including_K_1_slice))
    
        
        saddle_index_densities_excluding_K_1_slice = saddle_index_densities_including_K_1_slice[round.((saddle_index_densities_including_K_1_slice*normalization_factor)/neighbours_saddle_indices_shared_between) .!= 1.0]
        average_saddle_index_densities_slice_excluding_K_1_slice[i] = mean(saddle_index_densities_excluding_K_1_slice)

    end


    ## -- SWAP CALCULATION --

    # Get the unique E0 values
    E0_bin_values_swap = unique(energy_saddle_index_densities_data_matrix_swap[:,1])

    # For each E0 value, find all saddle index densities and average them
    average_saddle_index_densities_swap = zeros(Float64, length(E0_bin_values_swap))
    sem_saddle_index_densities_swap = zeros(Float64, length(E0_bin_values_swap))
    for (i, E0) in pairs(E0_bin_values_swap)
        # Get the saddle index densities for the current E0 swap
        saddle_index_densities_swap = [energy_saddle_index_densities_data_matrix_swap[j,2] for j in 1:size(energy_saddle_index_densities_data_matrix_swap,1) if energy_saddle_index_densities_data_matrix_swap[j,1] == E0] 
        average_saddle_index_densities_swap[i] = mean(saddle_index_densities_swap)
        sem_saddle_index_densities_swap[i] = std(saddle_index_densities_swap) / sqrt(length(saddle_index_densities_swap))
    end


    ### --- PLOTTING ---

    ## -- SWAP PLOTTING --

    # Plot the average saddle index densities against E0
    # But plot all values where the average saddle index = 0.0 in green
    zero_saddle_index_densities_swap_indices = [i for (i, saddle_index_density) in pairs(average_saddle_index_densities_swap) if saddle_index_density == 0.0]
    non_zero_saddle_index_densities_indices_swap = [i for i in 1:length(E0_bin_values_swap) if average_saddle_index_densities_swap[i] != 0.0]
    println(zero_saddle_index_densities_swap_indices)

    average_saddle_index_densities_graph = plot(
        legend=:topleft,
        xguidefontsize=12,  # X-axis label font size
        yguidefontsize=12,  # Y-axis label font size
        margin=1mm,         # Margin around the plot
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Average Saddle Index Density, "*L"\overline{\langle k \rangle}", 
        title = neighbour_order_to_measure_to>1 ? "Second Nearest Neighbours" : "",
        titlefontsize=12,
    )


    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube)), 
    average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap], 
    yerror=sem_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap],
    label="Swap-Move Cube",
    color=alex_blue,
    )

    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_swap[zero_saddle_index_densities_swap_indices]/-solved_configuration_energy(cube)), 
    zeros(Float64, length(zero_saddle_index_densities_swap_indices)), 
    yerror=sem_saddle_index_densities_swap[zero_saddle_index_densities_swap_indices],
    label="",
    color=alex_green,
    )

    if fitting 
        # Do exponential fit Ae^(bE) on 
        # average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap] against
        # -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))
        # and plot the fit
        fit = curve_fit((x, p) -> p[1]*exp.(p[2]*x), 
            -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube)), 
            average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap], 
            [1.0, 1.0])

        test_energies = range(minimum(-1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))), 
            maximum(-1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))), length=1000)

        plot!(average_saddle_index_densities_graph,
            test_energies, 
            fit.param[1]*exp.(fit.param[2]*test_energies),
            # label="$(round(fit.param[1], digits=2))e^{ $(round(fit.param[2], digits=2)) \\epsilon}",
            label="",
            color=:red,
            linestyle=:dash,
            linewidth=3)
            
        println("Exponential Fit Parameters For Swap-Move Cube: ", fit.param)


    end

    ## -- SLICE PLOTTING --

    # Plot the average saddle index densities against E0
    # But plot all values where the average saddle index = 0.0 in red
    # And plot all values where the average saddle index = 0.0 if you exclude K=1 saddles in orange
    zero_saddle_index_densities_including_K_1_indices = [i for (i, saddle_index_density_including_K_1) in pairs(average_saddle_index_densities_slice_including_K_1_slice) if saddle_index_density_including_K_1 == 0.0]
    zero_saddle_index_densities_excluding_K_1_indices = [i for (i, saddle_index_density_excluding_K_1) in pairs(average_saddle_index_densities_slice_excluding_K_1_slice) if saddle_index_density_excluding_K_1 == 0.0 && average_saddle_index_densities_slice_including_K_1_slice[i] != 0.0]
    non_zero_saddle_index_densities_indices_slice = [i for i in 1:length(E0_bin_values_slice) if average_saddle_index_densities_slice_including_K_1_slice[i] != 0.0 && average_saddle_index_densities_slice_excluding_K_1_slice[i] != 0.0]

    println("Slice zero saddle index densities including K=1 indices: ", zero_saddle_index_densities_including_K_1_indices)
    println("Slice zero saddle index densities excluding K=1 indices: ", zero_saddle_index_densities_excluding_K_1_indices)


    scatter!(average_saddle_index_densities_graph,
        -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube)), 
        average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice], 
        yerror=sem_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice],
        label="Slice-Rotation Cube",
        color=alex_orange,
    )

    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_slice[zero_saddle_index_densities_excluding_K_1_indices]/-solved_configuration_energy(cube)), 
    average_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_excluding_K_1_indices], 
    yerror=sem_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_excluding_K_1_indices],
    label="",
    color=alex_pink,
    )

    # Print largest E0 value for zero_saddle_index_densities_excluding_K_1_indices and then print E^*
    # println("Largest E0 value for zero saddle index densities excluding K=1 indices: ", -1.0 + maximum(E0_bin_values_slice[zero_saddle_index_densities_excluding_K_1_indices])/-solved_configuration_energy(cube))
    println("E^* = ", E_star)

    # Print lowest energy average saddle index density for each cube
    println("Lowest energy average saddle index density for swap-move cube: ", minimum(average_saddle_index_densities_swap))
    println("Lowest energy average saddle index density for slice-rotation cube: ", minimum(average_saddle_index_densities_slice_including_K_1_slice))

    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_slice[zero_saddle_index_densities_including_K_1_indices]/-solved_configuration_energy(cube)), 
    zeros(Float64, length(zero_saddle_index_densities_including_K_1_indices)), 
    yerror=sem_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_including_K_1_indices],
    label="",
    color=alex_red,
    )

    # Add E^* vertical line
    vline!(average_saddle_index_densities_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
    annotate!(average_saddle_index_densities_graph, [(E_star+0.025, ylims(average_saddle_index_densities_graph)[1]+0.66, Plots.text(L"\bar\epsilon^*", 12, :black))])

    # Add E_on vertical line
    vline!(average_saddle_index_densities_graph, [E_on], linecolor=alex_red, linestyle=:dash, linewidth=2, label="")
    annotate!(average_saddle_index_densities_graph, [(E_on+0.025, ylims(average_saddle_index_densities_graph)[1]+0.66, Plots.text(L"\bar\epsilon^{\rm on}", 12, :black))])

    if fitting
        # Do exponential fit Ae^(bE) on 
        # average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice] against
        # -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))
        # and plot the fit
        fit = curve_fit((x, p) -> p[1]*exp.(p[2]*x), 
            -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube)), 
            average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice], 
            [1.0, 1.0])

        test_energies = range(minimum(-1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))), 
            maximum(-1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))), length=1000)

        plot!(average_saddle_index_densities_graph,
            test_energies, 
            fit.param[1]*exp.(fit.param[2]*test_energies),
            # label="$(round(fit.param[1], digits=2))e^{ $(round(fit.param[2], digits=2)) \\epsilon}",
            label="",
            color=:black,
            linestyle=:dash,
            linewidth=3)
            
        println("Exponential Fit Parameters For Slice-Rotation Cube: ", fit.param)


    end


    #####


    ### --- LOG-LINEAR INSET GRAPH ---

    # Main inset data - swap data
    # Note: For log-scale error bars, we convert SEM to relative error by dividing by the mean
    # This is because d(ln(y))/dy = 1/y, so the error in ln(y) is approximately (error in y)/y
    scatter!(average_saddle_index_densities_graph, 
    -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]./-solved_configuration_energy(cube)),
    log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap]); 
    # yerror=sem_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap] ./ average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap],
    color=alex_blue, legend=false, inset=bbox(0.25,0.25,0.3,0.4), subplot=2,
    xlabel=L"\epsilon", ylabel=L"\ln\overline{\langle k \rangle}", yguidefontsize=12,xguidefontsize=12, ylims=(-11,0))

    # Main inset data - slice data
    scatter!(average_saddle_index_densities_graph, 
    -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube)),
    log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice]); 
    # yerror=sem_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice] ./ average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice],
    color=alex_orange, subplot=2)

    println("Maximum of log slice data: ", maximum(log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice])))
    println("Maximum of log swap data: ", maximum(log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap])))
    println("Minimum of log slice data: ", minimum(log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice])))
    println("Minimum of log swap data: ", minimum(log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap])))

    # Pink slice inset data
    # Same relative error conversion for log scale
    scatter!(average_saddle_index_densities_graph, 
        -1.0 .+ (E0_bin_values_slice[zero_saddle_index_densities_excluding_K_1_indices]/-solved_configuration_energy(cube)),
        log.(average_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_excluding_K_1_indices]); 
        # yerror=sem_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_excluding_K_1_indices] ./ average_saddle_index_densities_slice_including_K_1_slice[zero_saddle_index_densities_excluding_K_1_indices],
        color=alex_pink, subplot=2)
   

    # Add E^* vertical line
    vline!(average_saddle_index_densities_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="", subplot=2)
    # annotate!(average_saddle_index_densities_graph, [(E_star+0.025, ylims(average_saddle_index_densities_graph)[1]+0.66, Plots.text(L"E^*", 12, :black))])
       
    # Add E_on vertical line
    vline!(average_saddle_index_densities_graph, [E_on], linecolor=alex_red, linestyle=:dash, linewidth=2, label="", subplot=2)

    if fitting
        # Do linear fit Ax + B on 
        # log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap]) against
        # -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))
        # and plot the fit
        fit = curve_fit((x, p) -> p[1]*x .+ p[2], 
            -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube)), 
            log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap]), 
            [1.0, 1.0])

        test_energies = range(minimum(-1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))),
            maximum(-1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube))), length=1000)

        plot!(average_saddle_index_densities_graph,
            test_energies, 
            fit.param[1]*test_energies .+ fit.param[2],
            # label="$(round(fit.param[1], digits=2))x + $(round(fit.param[2], digits=2))",
            label="",
            color=:red,
            linestyle=:dash,
            linewidth=3, subplot=2)
        
        println("Linear Fit Parameters For Swap-Move Cube: ", fit.param)

        # Also plot on main graph as exponential
        # plot!(average_saddle_index_densities_graph,
        #     test_energies, 
        #     exp.(fit.param[1]*test_energies .+ fit.param[2]),
        #     # label="$(round(fit.param[1], digits=2))e^{ $(round(fit.param[2], digits=2)) \\epsilon}",
        #     label="",
        #     color=:red,
        #     linestyle=:dash,
        #     linewidth=3)

    end

    if fitting
        # Do linear fit Ax + B on 
        # log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice]) against
        # -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))
        # and plot the fit
        fit = curve_fit((x, p) -> p[1]*x .+ p[2], 
            -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube)), 
            log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice]), 
            [1.0, 1.0])

        test_energies = range(minimum(-1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))),
            maximum(-1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))), length=1000)

        plot!(average_saddle_index_densities_graph,
            test_energies, 
            fit.param[1]*test_energies .+ fit.param[2],
            # label="$(round(fit.param[1], digits=2))x + $(round(fit.param[2], digits=2))",
            label="",
            color=:black,
            linestyle=:dash,
            linewidth=3, subplot=2)
        
        println("Linear Fit Parameters For Slice-Rotation Cube: ", fit.param)

        # Also plot on main graph as exponential
        # plot!(average_saddle_index_densities_graph,
        #     test_energies, 
        #     exp.(fit.param[1]*test_energies .+ fit.param[2]),
        #     # label="$(round(fit.param[1], digits=2))e^{ $(round(fit.param[2], digits=2)) \\epsilon}",
        #     label="",
        #     color=:black,
        #     linestyle=:dash,
        #     linewidth=3)

    end

    ### --- SAVE AND DISPLAY GRAPH ---
    extra = fitting ? "_fitting" : ""

    savefig(average_saddle_index_densities_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_average_saddle_index_densities$(extra).png")
    savefig(average_saddle_index_densities_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_average_saddle_index_densities$(extra).pdf")
    display(average_saddle_index_densities_graph)
end


saddle_index_density_figure()

saddle_index_density_figure(fitting=true)