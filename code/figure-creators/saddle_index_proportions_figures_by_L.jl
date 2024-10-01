using DelimitedFiles
using Plots

# using StatsBase
using LaTeXStrings
using CSV
# using DataFrames

# using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

using LsqFit  
using Roots


include("../core/rubiks_cube.jl")

# function saddle_index_proportions_figure_by_L(L::Int64; neighbour_order_to_measure_to::Int64=1)

#     # E_star = -0.39015151515151514 
#     # E_star = -0.376098787878788
#     # E_star = -0.3859651515151515
#     # E_on = -0.23134348484848488
#  E_on -0.24010378787878783

#     simulation_name = "combined_disorder_average_connections_L=$(L)_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"

#     ### -- SET UP DEFAULT PARAMETERS --
#     header_line = readlines(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data",simulation_name))[1]
#     match_obj = match(r"L=(\d+)", header_line)
#     L = parse(Int, match_obj.captures[1])
#     cube = RubiksCube(L)

#     ### --- READ IN DATA ---
#     energy_saddle_index_densities_data_matrix = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data",simulation_name), ',', Float64, '\n', skipstart=3)


#     # Remove any rows with energy values larger than e10
#     energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .< 1e10, :]
#     # or with energy values >=0.0
#     energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(L)), :]


#     ### --- COLOURS ---
#     alex_red = RGB(227/255, 11/255, 92/255)
#     alex_pink = RGB(255/255, 105/255, 180/255)
#     alex_orange = RGB(255/255, 165/255, 0/255)
#     alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
#     alex_blue = RGB(100/255, 149/255, 237/255)


#     ### --- CREATE GRAPH ---
#     E0_bin_values = unique(energy_saddle_index_densities_data_matrix[:,1])

    
#     # Now gather the proportion of minima (saddle index K=0), K=1, K=2  and K>2 for each energy value
#     # and store in a new matrix against energy
#     k_saddle_proportions = zeros(Float64, length(E0_bin_values), 4)

#     Z = 6*(L-1)
#     normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
#     neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)

#     K_60 = 0

#     total = 0
#     for (i, E0) in pairs(E0_bin_values)
#         # Get the saddle indices for the current E0 slice
#         E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[j,2])/(neighbours_saddle_indices_shared_between)) for j in 1:size(energy_saddle_index_densities_data_matrix,1) if energy_saddle_index_densities_data_matrix[j,1] == E0] 
        
#         K_60 += sum(E0_saddle_indices .== 60)
#         total += length(E0_saddle_indices)


#         # Calculate the proportion of minima
#         k_saddle_proportions[i, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
#         k_saddle_proportions[i, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
#         k_saddle_proportions[i, 3] = sum(E0_saddle_indices .== 2) / length(E0_saddle_indices)
#         k_saddle_proportions[i, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
#     end

#     println("K=60 Maxima Saddles: ", K_60)
#     println("Total Saddles: ", total)

#     # Plot the k saddle proportions against E0
#     # Empty scatter just to put K=0 in correct legend and color position
#     k_saddle_proportions_graph = scatter([], [], 
#         label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)", 
#         color=alex_blue,
#         xlims=(-0.7, -0.05)
#     )

#     scatter!(k_saddle_proportions_graph,
#         -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
#     k_saddle_proportions[:, 4],
#     label="",
#     color=alex_green,
#     )

#     scatter!(k_saddle_proportions_graph,
#         -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
#         k_saddle_proportions[:, 2],
#         label="K=$(1+neighbour_order_to_measure_to-1)",
#         color=alex_pink,
#     )

#     scatter!(k_saddle_proportions_graph,
#         -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
#         k_saddle_proportions[:, 3],
#         label="K=$(2+neighbour_order_to_measure_to-1)",
#         color=alex_orange,
#     )

#     scatter!(k_saddle_proportions_graph,
#         -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
#         k_saddle_proportions[:, 1],
#         xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
#         ylabel="Saddle Index, "*L"K"*", Proportions",
#         # title="L=$L $connectivity Cube Saddle Index Proportions",
#         # titlefontsize=10,   # Title font size

#         # label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)",
#         label="",
#         xguidefontsize=12,   # X-axis label font size
#         yguidefontsize=12,   # Y-axis label font size
#         margin=0mm,          # Margin around the plot
#         legend=:left,
#         color=alex_blue,
#     )

#     # Empty scatter just to put K>2 in correct legend and color position
#     scatter!(k_saddle_proportions_graph, [], [], label="K≥$(2+neighbour_order_to_measure_to-1)", color=alex_green)


    
#     # TODO maybe add a vertical line at E^* and E^on later when know values for each L
#     # Add E^* vertical line
#     # if !isnothing(E_star)
#     #     vline!(k_saddle_proportions_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
#     #     annotate!(k_saddle_proportions_graph, [(E_star+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\epsilon^*", 12, :black))])
#     # end

#     # # Add E^on vertical line
#     # if !isnothing(E_on)
#     #     vline!(k_saddle_proportions_graph, [E_on], linecolor=:red, linestyle=:dash, linewidth=2, label="")
#     #     annotate!(k_saddle_proportions_graph, [(E_on+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\epsilon^{\rm on}", 12, :black))])
#     # end

#     # Save and display the graphs
#     # Remove string "_energy_saddle_index_densities" from simulation_name to get save_name
#     save_name = replace(simulation_name, "_energy_saddle_index_densities" => "")

#     savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/$(save_name)_k_saddle_proportions.pdf")
#     savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/$(save_name)_k_saddle_proportions.png")
#     display(k_saddle_proportions_graph)


# end











function all_saddle_index_proportions_figure_by_E1_E0_L_analysis(L_values::Array{Int64}; include_K1_in_K0::Bool=false)
    neighbour_order_to_measure_to::Int64=1

    addon = include_K1_in_K0 ? "_with_K1_in_K0" : ""

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    # Colors and symbols for different L values
    colors = [RGB(227/255, 11/255, 92/255), RGB(255/255, 165/255, 0/255), RGB(23/255,177/255,105/255), RGB(100/255, 149/255, 237/255), RGB(255/255, 105/255, 180/255)]



    ### --- READ IN DATA ---
    data_dict = Dict()

    for (i, L) in enumerate(L_values)
        simulation_name = "combined_disorder_average_connections_L=$(L)_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"
        file_path = joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data", simulation_name)
        
        cube = RubiksCube(L)
        
        energy_saddle_index_densities_data_matrix = readdlm(file_path, ',', Float64, '\n', skipstart=3)
        
        # Remove rows with energy values larger than e10 or >=0.0
        energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[
            (energy_saddle_index_densities_data_matrix[:,1] .< 1e10) .& 
            (energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(L))), :]
        
        E0_bin_values = unique(energy_saddle_index_densities_data_matrix[:,1])
        
        Z = 6*(L-1)
        normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
        neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)
        
        k_saddle_proportions = zeros(Float64, length(E0_bin_values), 4)
        
        for (j, E0) in pairs(E0_bin_values)
            E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[k,2])/(neighbours_saddle_indices_shared_between)) 
                                 for k in 1:size(energy_saddle_index_densities_data_matrix,1) 
                                 if energy_saddle_index_densities_data_matrix[k,1] == E0]
            
            k_saddle_proportions[j, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
            k_saddle_proportions[j, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
            
            if include_K1_in_K0
                k_saddle_proportions[j,1] = k_saddle_proportions[j,1] + k_saddle_proportions[j,2]
            end
            
            k_saddle_proportions[j, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
        end
        
        data_dict[L] = (E0_bin_values, k_saddle_proportions)
    end

        # Dictionary for epsilon_star values
        epsilon_star = Dict(
            3 => -0.49,
            5 => -0.49,
            7 => -0.45,
            9 => -0.41,
            11 => -0.39
        )
    
        # Dictionary for epsilon_on values
        epsilon_on = Dict(
            3 => -0.40,
            # 5 => -0.31,
            5 => -0.24,
            7 => -0.30,
            9 => -0.28,
            11 => -0.23
        )









    ### --- CREATE INDIVIDUAL GRAPHS FOR EACH L ---
    intersection_points_dict = Dict()
    epsilon_times_dict = Dict()
    epsilon_otims_dict = Dict()

    cutoff_saddle_index_proportion = 1e-2


    for (i, L) in enumerate(L_values)
        individual_L_plot = scatter(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
            ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
            legend=:left,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05)
        )

        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        scatter!(individual_L_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 1],
            # label="K=0, L=$L",
            label="",
            color=colors[i],
            marker=:circle,
            markersize=4
        )

        scatter!(individual_L_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            # label="K≥2, L=$L",
            label="",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )

        # Also plot 1 - K2
        # scatter!(individual_L_plot,
        #     -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
        #     1 .- k_saddle_proportions[:, 4],
        #     # label="K≥2, L=$L",
        #     label="",
        #     color=colors[i],
        #     marker=:square,
        #     markersize=4,
        # )

        # Add E^* vertical line
        if haskey(epsilon_star, L)
            vline!(individual_L_plot, [epsilon_star[L]], linecolor=:green, linestyle=:dash, linewidth=2, label="")
            annotate!(individual_L_plot, [(epsilon_star[L]+0.023, ylims(individual_L_plot)[1]+0.58, Plots.text(L"\epsilon^*", 12, :black))])
        end

        # Add E^on vertical line
        if haskey(epsilon_on, L)
            vline!(individual_L_plot, [epsilon_on[L]], linecolor=:red, linestyle=:dash, linewidth=2, label="")
            annotate!(individual_L_plot, [(epsilon_on[L]+0.023, ylims(individual_L_plot)[1]+0.58, Plots.text(L"\epsilon^{\rm on}", 12, :black))])
        end
    








        # Run a LsQ fit on the K=0 and K>=2 data
        K0_functional_form(x, p) = (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^(1/p[3])
        K2_functional_form(x, p) = 1 .- K0_functional_form(x, p)
        p0 = [30.0, -0.7, 0.2]  # Initial guess for [k, x0, v]


        # Fit the K=0 data
        valid_rows_K0 = .!isnan.(k_saddle_proportions[:, 1]) .& .!isinf.(k_saddle_proportions[:, 1]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K0 = E0_bin_values[valid_rows_K0]
        k_saddle_proportions_filtered_K0 = k_saddle_proportions[valid_rows_K0, 1]
        
        # Fit the K>=2 data
        valid_rows_K2 = .!isnan.(k_saddle_proportions[:, 4]) .& .!isinf.(k_saddle_proportions[:, 4]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K2 = E0_bin_values[valid_rows_K2]
        k_saddle_proportions_filtered_K2 = k_saddle_proportions[valid_rows_K2, 4]
        
        # Transform E0_bin_values to energy density
        energy_density_K0 = -1.0 .+ (E0_bin_values_filtered_K0./-solved_configuration_energy(RubiksCube(L)))
        energy_density_K2 = -1.0 .+ (E0_bin_values_filtered_K2./-solved_configuration_energy(RubiksCube(L)))
        

        fit_K0 = curve_fit(K0_functional_form, energy_density_K0, k_saddle_proportions_filtered_K0, p0, show_trace=false, maxIter=200)
        fit_K2 = curve_fit(K2_functional_form, energy_density_K2, k_saddle_proportions_filtered_K2, p0, show_trace=false, maxIter=200)
        

        # Plot fitted curves
        fit_values = LinRange(-1.0, -0.05, 100)
        plot!(
            individual_L_plot,
            fit_values,
            K0_functional_form(fit_values, fit_K0.param),
            label="",
            color=colors[i],
            linestyle=:dash,
        )

        plot!(
            individual_L_plot,
            fit_values,
            K2_functional_form(fit_values, fit_K2.param),
            label="",
            color=colors[i],
            linestyle=:dash,
        )
        

        # # Find the intersection point
        intersection_difference(x) = K0_functional_form(x, fit_K0.param) - K2_functional_form(x, fit_K2.param)
        intersection_point = find_zero(intersection_difference, (-2.0, 0.0), Bisection())
        intersection_point_saddle_index_proportion = K0_functional_form(intersection_point, fit_K0.param)
        intersection_points_dict[L] = (intersection_point, intersection_point_saddle_index_proportion)
        vline!(individual_L_plot, [intersection_point], linecolor=:black, linestyle=:dash, linewidth=2, label="Intersection Point")
        
        println("K0 Fit: ", fit_K0.param)
        println("K2 Fit: ", fit_K2.param)
        println("L=$L Intersection Point: ", intersection_point, " Saddle Index Proportion: ", intersection_point_saddle_index_proportion)


        # Find epsilon_times which is the energy density where K>=2 drops below cutoff_saddle_index_proportion
        # using the K>=2 fit
        epsilon_times_difference(x) = K2_functional_form(x, fit_K2.param) - cutoff_saddle_index_proportion
        epsilon_times = find_zero(epsilon_times_difference, (-2.0, 0.0), Bisection())
        epsilon_times_dict[L] = epsilon_times
        vline!(individual_L_plot, [epsilon_times], linecolor=:pink, linestyle=:dash, linewidth=2, label=L"\bar \epsilon_{\times}")
        println("L=$L epsilon_\times: ", epsilon_times)

        # Find epsilon_otims which is the energy density where K=0 rises above cutoff_saddle_index_proportion
        # using the K=0 fit
        epsilon_otimes_difference(x) = K0_functional_form(x, fit_K0.param) - cutoff_saddle_index_proportion
        epsilon_otims = find_zero(epsilon_otimes_difference, (-2.0, 0.0), Bisection())
        epsilon_otims_dict[L] = epsilon_otims
        vline!(individual_L_plot, [epsilon_otims], linecolor=:orange, linestyle=:dash, linewidth=2, label=L"\bar \epsilon_{\otimes}")
        println("L=$L epsilon_otimes: ", epsilon_otims)

        savefig(individual_L_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_k_saddle_proportions.pdf")
        savefig(individual_L_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_k_saddle_proportions.png")
        display(individual_L_plot)
    end


    println("Intersection points dict: ", intersection_points_dict)
    println("Epsilon times dict: ", epsilon_times_dict)
    println("Epsilon otims dict: ", epsilon_otims_dict)








end



function combined_saddle_index_proportions_figure_by_L(L_values::Array{Int64}; include_K1_in_K0::Bool=false)
    neighbour_order_to_measure_to::Int64=1

    addon = include_K1_in_K0 ? "_with_K1_in_K0" : ""

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    # Colors and symbols for different L values
    colors = [RGB(227/255, 11/255, 92/255), RGB(255/255, 165/255, 0/255), RGB(23/255,177/255,105/255), RGB(100/255, 149/255, 237/255), RGB(255/255, 105/255, 180/255)]



    ### --- READ IN DATA ---
    data_dict = Dict()

    epsilon_log_dict = Dict()

    for (i, L) in enumerate(L_values)
        simulation_name = "combined_disorder_average_connections_L=$(L)_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"
        file_path = joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data", simulation_name)
        
        cube = RubiksCube(L)
        
        energy_saddle_index_densities_data_matrix = readdlm(file_path, ',', Float64, '\n', skipstart=3)
        
        # Remove rows with energy values larger than e10 or >=0.0
        energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[
            (energy_saddle_index_densities_data_matrix[:,1] .< 1e10) .& 
            (energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(L))), :]
        
        E0_bin_values = unique(energy_saddle_index_densities_data_matrix[:,1])
        
        Z = 6*(L-1)
        normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
        neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)
        
        k_saddle_proportions = zeros(Float64, length(E0_bin_values), 4)
        
        for (j, E0) in pairs(E0_bin_values)
            E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[k,2])/(neighbours_saddle_indices_shared_between)) 
                                 for k in 1:size(energy_saddle_index_densities_data_matrix,1) 
                                 if energy_saddle_index_densities_data_matrix[k,1] == E0]
            
            k_saddle_proportions[j, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
            k_saddle_proportions[j, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
            
            if include_K1_in_K0
                k_saddle_proportions[j,1] = k_saddle_proportions[j,1] + k_saddle_proportions[j,2]
            end
            
            k_saddle_proportions[j, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
        end
        
        data_dict[L] = (E0_bin_values, k_saddle_proportions)
    end







    ### --- CREATE COMBINED GRAPH ---
    combined_plot = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
        # legend=:topright,
        legend=(0.9,0.8),
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05)
    )


    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        scatter!(combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 1],
            # label="K=0, L=$L",
            label="",
            color=colors[i],
            marker=:circle,
            markersize=4
        )

        scatter!(combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            # label="K≥2, L=$L",
            label="",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )
    end






    # Run a LsQ fit on the K=0 and K>=2 data
    K0_functional_form(x, p) = (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^(1/p[3])
    K2_functional_form(x, p) = 1 .- K0_functional_form(x, p)
    p0 = [30.0, -0.7, 0.2]  # Initial guess for [k, x0, v]


    # Fit the data
    fit_data_dict = Dict()

    L_values = [11, 9, 7, 5, 3]
    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Fit the K=0 data
        valid_rows_K0 = .!isnan.(k_saddle_proportions[:, 1]) .& .!isinf.(k_saddle_proportions[:, 1]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K0 = E0_bin_values[valid_rows_K0]
        k_saddle_proportions_filtered_K0 = k_saddle_proportions[valid_rows_K0, 1]
        
        # Fit the K>=2 data
        valid_rows_K2 = .!isnan.(k_saddle_proportions[:, 4]) .& .!isinf.(k_saddle_proportions[:, 4]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K2 = E0_bin_values[valid_rows_K2]
        k_saddle_proportions_filtered_K2 = k_saddle_proportions[valid_rows_K2, 4]
        
        # Transform E0_bin_values to energy density
        energy_density_K0 = -1.0 .+ (E0_bin_values_filtered_K0./-solved_configuration_energy(RubiksCube(L)))
        energy_density_K2 = -1.0 .+ (E0_bin_values_filtered_K2./-solved_configuration_energy(RubiksCube(L)))
        

        fit_K0 = curve_fit(K0_functional_form, energy_density_K0, k_saddle_proportions_filtered_K0, p0, show_trace=true, maxIter=200)
        fit_K2 = curve_fit(K2_functional_form, energy_density_K2, k_saddle_proportions_filtered_K2, p0, show_trace=true, maxIter=200)
        
        fit_data_dict[L] = (fit_K0, fit_K2)

        println("L=$L")
        println("K=0 Fit: ", fit_K0.param)
        println("K>=2 Fit: ", fit_K2.param)
    end




    # Plot individual fits on graph and also extract intersection points
    intersection_energy_density = []
    intersection_saddle_index_proportion = []

    collapsed_plot = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Collapsed "*L"K \geq 2"*" Proportions",
        legend=:left,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05)
    )

    # Initialize a new plot where x-axis is alpha*(epsilon - epsilon_log)
    energy_density_transformed_plot = scatter(
        xlabel=L"\alpha (\epsilon - \epsilon_{\log})",
        ylabel=L"K \geq 2"*" Proportions",
        legend=:left,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm
    )

    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        fit_K0, fit_K2 = fit_data_dict[L]
    
        # Plot fitted curves
        fit_values = LinRange(-1.0, -0.05, 100)
        # plot!(
        #     combined_plot,
        #     fit_values,
        #     K0_functional_form(fit_values, fit_K0.param),
        #     label="",
        #     color=colors[i],
        #     linestyle=:dash,
        # )

        plot!(
            combined_plot,
            fit_values,
            K2_functional_form(fit_values, fit_K2.param),
            label="",
            color=colors[i],
            linestyle=:dash,
        )


        # Display temporary graph with individual L data and fit
        individual_L_plot = scatter(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
            ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
            legend=:left,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05)
        )
        scatter!(individual_L_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            label="K≥2, L=$L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )
        fit_values = LinRange(-1.0, -0.05, 100)
        plot!(
            individual_L_plot,
            fit_values,
            K2_functional_form(fit_values, fit_K2.param),
            label="",
            color=colors[i],
            linestyle=:dash,
        )
        display(individual_L_plot)




        # If \gamma = 1/fit_K2.param[3] 
        # If \epsilon_log = fit_K2.param[2]
        # If \alpha = fit_K2.param[1]
        alpha = fit_K2.param[1]
        epsilon_log = fit_K2.param[2]
        gamma = 1/fit_K2.param[3]



        # COLLAPSED PLOT
        # Plot collapsed data for K>=2
        # f(\epsilon) = (1/alpha)* ln((1-pK2)^(1/gamma) - 1) + \epsilon_log
        # and it should all collapse to f(\epsilon) = \epsilon
        
        # Compute the argument of the log function
        log_args = (1 .- k_saddle_proportions[:, 4]) .^ (-(1/gamma)) .- 1

        # Identify positive arguments
        positive_indices = log_args .> 0

        # Initialize the collapsed values array
        collapsed_K2_values = similar(log_args)

        # Apply the formula only to positive arguments
        collapsed_K2_values[positive_indices] = (1 / alpha) * log.(log_args[positive_indices]) .+ epsilon_log


        # Handle negative or zero arguments (e.g., set to NaN or skip)
        collapsed_K2_values[.!positive_indices] .= NaN  # or handle as needed
        println("Not positive indices: ", .!positive_indices)
        
        
        scatter!(collapsed_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            collapsed_K2_values,
            label="",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )

        # Also plot black dashed y=x line 
        plot!(
            collapsed_plot,
            LinRange(-1.0, -0.05, 100),
            LinRange(-1.0, -0.05, 100),
            label="",
            color=:black,
            linestyle=:dash,
        )



        # TRANSFORMED ENERGY DENSITY PLOT
        # Calculate the energy density for K>=2
        energy_density_K2 = -1.0 .+ (E0_bin_values ./ -solved_configuration_energy(RubiksCube(L)))

        # Transform the energy density: alpha*(epsilon - epsilon_log)
        transformed_energy_density = alpha * (energy_density_K2 .- epsilon_log)

        # Plot the transformed energy density against K>=2 proportions
        scatter!(
            energy_density_transformed_plot,
            transformed_energy_density[abs.(transformed_energy_density) .< 1e3],
            k_saddle_proportions[abs.(transformed_energy_density) .< 1e3, 4],  # Proportion of K >= 2
            label="L = $L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )


   end

   display(collapsed_plot)
   savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_collapsed$(addon).pdf")
   savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_collapsed$(addon).png")


    display(energy_density_transformed_plot)
    savefig(energy_density_transformed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_transformed_energy_density$(addon).pdf")
    savefig(energy_density_transformed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_transformed_energy_density$(addon).png")
















    L_values = [11, 9, 7, 5]

    ### --- BETA AND GAMMA FIT PARAMETERS GRAPH ---
    plot = scatter(
        # xlabel=L"L",
        xlabel=L"1/L",
        ylabel=L"\beta, \ \gamma"*" Fit Parameter Values",
        legend=:right,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
    )

    K0_fit_values = [fit_data_dict[L][1].param for L in L_values]
    K2_fit_values = [fit_data_dict[L][2].param for L in L_values]


    # scatter!(plot, [L for L in L_values], [x[3] for x in K2_fit_values], label="K≥2, "*L"\gamma", color=:red, marker=:diamond, markersize=4)
    # scatter!(plot, [L for L in L_values], [x[3] for x in K0_fit_values], label="K=0, "*L"\gamma", color=:red, marker=:square, markersize=4)
    scatter!(plot, [1/L for L in L_values], [x[3] for x in K2_fit_values], label="K≥2, "*L"\gamma", color=:red, marker=:diamond, markersize=4)
    # scatter!(plot, [1/L for L in L_values], [x[3] for x in K0_fit_values], label="K=0, "*L"\gamma", color=:red, marker=:square, markersize=4)
   


    # scatter!(plot, [L for L in L_values], [x[2] for x in K0_fit_values], label="K=0, "*L"\beta", color=:blue, marker=:square, markersize=4)
    # scatter!(plot, [L for L in L_values], [x[2] for x in K2_fit_values], label="K≥2, "*L"\beta", color=:blue, marker=:diamond, markersize=4)
    # scatter!(plot, [1/L for L in L_values], [x[2] for x in K0_fit_values], label="K=0, "*L"\beta", color=:blue, marker=:square, markersize=4)
    scatter!(plot, [1/L for L in L_values], [x[2] for x in K2_fit_values], label="K≥2, "*L"\beta", color=:blue, marker=:diamond, markersize=4)



    # fit_beta_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_values, [x[2] for x in K0_fit_values], [0.0, 0.0])
    # m, c = fit_beta_K0.param
    # confidence_intervals = confidence_interval(fit_beta_K0, 0.1)
    # c_ci = confidence_intervals[2]
    # c_err = abs(c - c_ci[1])
    # fitting_1_L_values = LinRange(0, 0.2, 100)
    # plot!(plot, fitting_1_L_values, m .+ c.*fitting_1_L_values, label="$(round(m, digits=2)) + $(round(c, digits=2))L, K=0", color=:black, linestyle=:dash)


    display(plot)
    savefig(plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_fit_parameters_beta_gamma$(addon).pdf")
    savefig(plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_fit_parameters_beta_gamma$(addon).png")





    ### --- ALPHA PARAMETERS GRAPH ---
    plot = scatter(
        xlabel=L"L",
        ylabel=L"\alpha"*" Fit Parameter Value",
        legend=:topleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
    )

    scatter!(plot, [L for L in L_values], [x[1] for x in K0_fit_values], label="K=0", color=:green, marker=:square, markersize=4)
    scatter!(plot, [L for L in L_values], [x[1] for x in K2_fit_values], label="K=2", color=:green, marker=:diamond, markersize=4)

    # Add linear fits for alpha values
    fit_alpha_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_values, [x[1] for x in K0_fit_values], [0.0, 0.0])
    fit_alpha_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_values, [x[1] for x in K2_fit_values], [0.0, 0.0])

    plot!(plot, L_values, fit_alpha_K0.param[1] .+ fit_alpha_K0.param[2].*L_values, label="$(round(fit_alpha_K0.param[1], digits=2)) + $(round(fit_alpha_K0.param[2], digits=2))L, K=0", color=:black, linestyle=:dash)
    plot!(plot, L_values, fit_alpha_K2.param[1] .+ fit_alpha_K2.param[2].*L_values, label="$(round(fit_alpha_K2.param[1], digits=2)) + $(round(fit_alpha_K2.param[2], digits=2))L, K≥2", color=:black, linestyle=:dash)
    
    display(plot)
    savefig(plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_fit_parameters_alpha$(addon).pdf")
    savefig(plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_fit_parameters_alpha$(addon).png")





















    ### --- K0 AND K2 FLOW PLOTS ---


    # Choose 10 values between energy density -0.5 and -0.1, and use the fits to plot their value as 1/L changes
    energy_density_values = vcat(LinRange(-0.26-0.2, -0.265, 15), LinRange(-0.25995, -0.23, 3), 0.15)
    L_values = [11,9,7,5]

    # Find the value in energy_density closest to -0.26
    closest_value = energy_density_values[argmin(abs.(energy_density_values .+ 0.26))]


    ## K0 PLOT
    K0_plot = scatter(
        xlabel=L"1/L",
        ylabel=L"K=0"*" Minima Proportions",
        legend=:topleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
    )



    plot!(K0_plot, [], [], color=alex_blue, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon < -0.26")
    for (i, energy_density) in pairs(energy_density_values)
        if energy_density == closest_value
            linestyle = :dash
            color = :black
            label = L"$\epsilon = -0.26$"
        else
            linestyle = :solid
            # color = colors[mod(i,length(colors))+1]
            color = energy_density > closest_value ? alex_red : alex_blue
            label=""
        end
        plot!(K0_plot, [1/L for L in L_values], [K0_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], color=color, linewidth=3, legend=:right, linestyle=linestyle, label=label)


        # Also fit the data and plot the fit
        
        if energy_density <= closest_value
            one_minus_fitting_function(x,p) =  1 .- p[1].*exp.(p[2] .* x)
            fitting_1_L_values = LinRange(0, 0.2, 100)
            fit = curve_fit(one_minus_fitting_function, L_values, [K0_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], [0.0, 0.0], show_trace=false, maxIter=200)
            
            # If fit has value that exceeds 1 or -1 then do nothing, otherwise plot
            if all(one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param) .<= 1) && all(one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param) .>= -1)
                plot!(K0_plot, fitting_1_L_values, one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param), color=alex_blue, linestyle=:dash, label="")
            end

        else
            fitting_function(x,p) =  p[1].*exp.(p[2] .* x)
            fitting_1_L_values = LinRange(0, 0.2, 100)
            fit = curve_fit(fitting_function, L_values, [K0_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], [0.0, 0.0], show_trace=false, maxIter=200)

            # If fit has value that exceeds 1 or -1 then do nothing, otherwise plot
            if all(fitting_function(1 ./ fitting_1_L_values, fit.param) .<= 1) && all(fitting_function(1 ./ fitting_1_L_values, fit.param) .>= -1)
                plot!(K0_plot, fitting_1_L_values, fitting_function(1 ./ fitting_1_L_values, fit.param), color=alex_red, linestyle=:dash, label="")
            end
        end
        
    end
    plot!(K0_plot, [], [], color=alex_red, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon > -0.26")

    display(K0_plot)
    savefig(K0_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_K0_flow$(addon).pdf")
    savefig(K0_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_K0_flow$(addon).png")





    ## K2 PLOT
    K2_plot = scatter(
        xlabel=L"1/L",
        ylabel=L"K\geq2"*" Saddle Proportions",
        legend=:topleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
    )

    plot!(K2_plot, [], [], color=alex_red, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon > -0.26")
    for (i, energy_density) in pairs(energy_density_values)
        if energy_density == closest_value
            linestyle = :dash
            color = :black
            label=L"$\epsilon = -0.26$"
        else
            linestyle = :solid
            color = energy_density > closest_value ? alex_red : alex_blue
            label=""
        end
        plot!(K2_plot, [1/L for L in L_values], [K2_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], color=color, linewidth=3, legend=:right, linestyle=linestyle, label=label)

        # Also fit the data and plot the fit
        
        if energy_density <= closest_value
            fitting_function(x,p) =  p[1].*exp.(p[2] .* x)
            fitting_1_L_values = LinRange(0, 0.2, 100)
            fit = curve_fit(fitting_function, L_values, [K2_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], [0.0, 0.0], show_trace=false, maxIter=200)
            
            # If fit has value that exceeds 1 or -1 then do nothing, otherwise plot
            if all(fitting_function(1 ./ fitting_1_L_values, fit.param) .<= 1) && all(fitting_function(1 ./ fitting_1_L_values, fit.param) .>= -1)
                plot!(K2_plot, fitting_1_L_values, fitting_function(1 ./ fitting_1_L_values, fit.param), color=alex_blue, linestyle=:dash, label="")
            end
        else
            one_minus_fitting_function(x,p) =  1 .- p[1].*exp.(p[2] .* x)
            fitting_1_L_values = LinRange(0, 0.2, 100)
            fit = curve_fit(one_minus_fitting_function, L_values, [K2_functional_form(energy_density, fit_data_dict[L][1].param) for (i, L) in enumerate(L_values)], [0.0, 0.0], show_trace=false, maxIter=200)

            # If fit has value that exceeds 1 or -1 then do nothing, otherwise plot
            if all(one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param) .<= 1) && all(one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param) .>= -1)
                plot!(K2_plot, fitting_1_L_values, one_minus_fitting_function(1 ./ fitting_1_L_values, fit.param), color=alex_red, linestyle=:dash, label="")
            end
        end
    end
    plot!(K2_plot, [], [], color=alex_blue, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon < -0.26")

    display(K2_plot)
    savefig(K2_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_K2_flow$(addon).pdf")
    savefig(K2_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_K2_flow$(addon).png")





    # Add alpha parameter graph to inset of combined plot
    # scatter!(combined_plot, [L for L in L_values], [x[1] for x in K0_fit_values], label="K=0", color=:green, marker=:circle, markersize=4, inset=bbox(0.24,0.225,0.2,0.2), subplot=2, ylabel=L"\alpha", xlabel=L"L", legend=false)
    scatter!(combined_plot, [L for L in L_values], [x[1] for x in K2_fit_values], label="K=2", color=:green, marker=:diamond, markersize=4, inset=bbox(0.25,0.225,0.2,0.2), subplot=2, ylabel=L"\alpha", xlabel=L"L", legend=false, ylims=(20,42), xlims=(4.5,11.5))
    # plot!(combined_plot, L_values, fit_alpha_K0.param[1] .+ fit_alpha_K0.param[2].*L_values, label="$(round(fit_alpha_K0.param[1], digits=2)) + $(round(fit_alpha_K0.param[2], digits=2))L, K=0", color=:black, linestyle=:dash, subplot=2)
    plot!(combined_plot, L_values, fit_alpha_K2.param[1] .+ fit_alpha_K2.param[2].*L_values, label="", color=:black, linestyle=:dash, subplot=2)
    
    # Add beta and gamma parameter graph to inset of combined plot
    # scatter!(combined_plot, [L for L in L_values], [x[3] for x in K2_fit_values], label="K≥2, "*L"\gamma", color=:red, marker=:diamond, markersize=4, inset=bbox(0.24,0.515,0.2,0.2), subplot=3, ylabel=L"\beta, \, \gamma", xlabel=L"L", legend=false, ylims=(-0.25,0.25))
    # scatter!(combined_plot, [L for L in L_values], [x[3] for x in K0_fit_values], label="K=0, "*L"\gamma", color=:red, marker=:circle, markersize=4, subplot=3)
    
    # scatter!(combined_plot, [L for L in L_values], [1/x[3] for x in K2_fit_values], label="K≥2, "*L"\gamma", color=:red, marker=:diamond, markersize=4, inset=bbox(0.24,0.515,0.2,0.2), subplot=3, ylabel=L"\beta, \, \gamma", xlabel=L"L", legend=false)
    # scatter!(combined_plot, [L for L in L_values], [1/x[3] for x in K0_fit_values], label="K=0, "*L"\gamma", color=:red, marker=:circle, markersize=4, subplot=3)
    

    scatter!(combined_plot, [1/L for L in L_values], [1/x[3] for x in K2_fit_values], label="K≥2, "*L"\gamma", color=:red, marker=:diamond, markersize=4, inset=bbox(0.75,0.515,0.2,0.2), subplot=3, ylabel=L"\gamma", xlabel=L"1/L", legend=false, ylims=(5.0,8.0), xlims=(0.075,0.21), xticks=[0.1,0.15,0.2]) 




    
    # scatter!(combined_plot, [L for L in L_values], [x[2] for x in K0_fit_values], label="K=0, "*L"\beta", color=:blue, marker=:circle, markersize=4, subplot=3)
    # scatter!(combined_plot, [L for L in L_values], [x[2] for x in K2_fit_values], label="K≥2, "*L"\beta", color=:blue, marker=:diamond, markersize=4, subplot=3)

    scatter!(combined_plot, [1/L for L in L_values], [x[2] for x in K2_fit_values], label="", color=:blue, marker=:diamond, markersize=4, subplot=4, ylabel=L"\bar \bar\epsilon_{\rm log}", xlabel=L"1/L", inset=bbox(0.255,0.515,0.2,0.2), ylims=(-0.245,-0.20), xlims=(0.0,0.22), xticks=[0,0.1,0.2])
    # Add linear fits for beta values
    fit_beta_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, 1 ./ L_values, [x[2] for x in K2_fit_values], [0.0, 0.0])
    c, m  = fit_beta_K2.param
    confidence_intervals = confidence_interval(fit_beta_K2, 0.1)
    c_ci = confidence_intervals[2]
    c_error = c - c_ci[1]
    fitting_1_L_values = LinRange(0, 0.2, 100)
    plot!(combined_plot, fitting_1_L_values, c .+ m.*fitting_1_L_values, label="", color=:blue, linestyle=:dash, subplot=4)
    # scatter!(combined_plot, [0.0], [c], yerror=(c_error), subplot=3, color=:blue, markerstrokecolor=:blue, markersize=1, label="", errorbar_color=:blue)
    println("K=2 Beta Fit: ", fit_beta_K2.param)
    println("L values: ", L_values)
    println("\bar \epsilon_{log} values: ", [x[2] for x in K2_fit_values])




    # Build combined plot legend
    scatter!(combined_plot, [], [], label="K=0", color=:white, marker=:circle, markersize=4)
    scatter!(combined_plot, [], [], label="K≥2", color=:white, marker=:diamond, markersize=4)
    scatter!(combined_plot, [], [], label="L=11", color=colors[1], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=9", color=colors[2], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=7", color=colors[3], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=5", color=colors[4], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=3", color=colors[5], marker=:square, markersize=4)

    # Save and display the combined graph
    save_name = "combined_saddle_index_proportions_multiple_L"
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)$(addon).pdf")
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)$(addon).png")
    display(combined_plot)








end



# all_saddle_index_proportions_figure_by_E1_E0_L_analysis([11,9,7,5,3])


# combined_saddle_index_proportions_figure_by_L([11,9,7,5,3]; include_K1_in_K0=true)
combined_saddle_index_proportions_figure_by_L([11,9,7,5,3])






























