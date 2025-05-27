import Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

include("../core/rubiks_cube.jl")

using DelimitedFiles
using Plots
using LaTeXStrings
using CSV
using Statistics 
using Plots.PlotMeasures
using LsqFit  
using Roots
using Dierckx  # Add this for spline interpolation
using SpecialFunctions  # Add this for erf function



# Define logistic function with fixed epsilon_log
function K2_fixed_epsilon_log(x, p, epsilon_log_fixed)
    # p[1] = alpha, p[2] = gamma
    return 1 .- (1 ./ (1 .+ exp.(p[1] .* (x .- epsilon_log_fixed)))).^(p[2])
end

# Define logistic function with fixed epsilon_log
function K2_fixed_epsilon_log_and_gamma(x, p, epsilon_log_fixed, gamma_fixed)
    # p[1] = alpha, p[2] = gamma
    return 1 .- (1 ./ (1 .+ exp.(p[1] .* (x .- epsilon_log_fixed)))).^(gamma_fixed)
end












function all_saddle_index_proportions_figure_by_E1_E0_L_analysis(L_values::Array{Int64}; include_K1_in_K0::Bool=false)
    neighbour_order_to_measure_to::Int64=1

    addon = include_K1_in_K0 ? "_with_K1_in_K0" : ""

    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255) # TODO add alt blue

    # Colors and symbols for different L values
    colors = [RGB(227/255, 11/255, 92/255), RGB(255/255, 165/255, 0/255), RGB(23/255,177/255,105/255), RGB(100/255, 149/255, 237/255), RGB(255/255, 105/255, 180/255), RGB(159/255, 226/255, 191/255), RGB(128/255, 0/255, 128/255)]



    ### --- READ IN DATA ---
    data_dict = Dict()

    for (i, L) in enumerate(L_values)
        println("Processing L=$L")
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
        K0_functional_form(x, p) = (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^(p[3])  # p[3] is now gamma directly
        K2_functional_form(x, p) = 1 .- K0_functional_form(x, p)
        p0 = [30.0, -0.7, 5.0]  # Initial guess for [alpha, epsilon_log, gamma]


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
            color=:black,
            linestyle=:dash,
            linewidth=2
        )

        plot!(
            individual_L_plot,
            fit_values,
            K2_functional_form(fit_values, fit_K2.param),
            label="",
            color=:black,
            linestyle=:dash,
            linewidth=2
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
    Plots.default(dpi = 600)
    
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)


    # Colors and symbols for different L values
    colors = [alex_red, alex_orange, alex_green, alex_blue, alex_pink, RGB(159/255, 226/255, 191/255), RGB(128/255, 0/255, 128/255)]


    ##### READ IN DATA #####
    data_dict = Dict()

    epsilon_log_dict = Dict()

    for (i, L) in enumerate(L_values)
        println("Processing L=$L")
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

    ##### #####












    ##### CREATE COMBINED SADDLE INDEX PROPORTIONS GRAPH SCATTER #####
    combined_plot = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
        # legend=:topright,
        legend=(0.92,0.85),
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


    ##### #####












    ##### CREATE COMBINED JUST K=1 SADDLE INDEX PROPORTIONS GRAPH #####
    combined_plot_K1 = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=1"*" Proportions",
        # legend=:topright,
        legend=(0.9,0.8),
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05)
    )

    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        scatter!(combined_plot_K1,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 2],
            label="L=$L",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
    end

    savefig(combined_plot_K1, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K=1_saddle_proportions.pdf")
    savefig(combined_plot_K1, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K=1_saddle_proportions.png")
    display(combined_plot_K1)


    ##### #####








    ##### GENERALISED LOGISTIC FUNCTION K0 AND K2 FIT CALCULATION #####


    # Run a LsQ fit on the K=0 and K>=2 data
    K0_functional_form(x, p) = (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^(p[3])  # p[3] is now gamma directly
    K2_functional_form(x, p) = 1 .- K0_functional_form(x, p)
    p0 = [30.0, -0.7, 5.0]  # Initial guess for [alpha, epsilon_log, gamma]


    # Fit the data
    fit_data_dict = Dict()

    L_values = [15, 13, 11, 9, 7, 5, 3]
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

    ##### #####











    ##### GENERALISED LOGISTIC FUNCTION COLLAPSED PLOTS #####


    collapsed_plot = scatter(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Collapsed "*L"K \geq 2"*" Proportions",
        legend=:bottomright,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
        xlims=(-0.7, -0.05)
    )

    # Initialize a new plot where x-axis is alpha(L)*(epsilon - epsilon_log)
    energy_density_transformed_plot = scatter(
        xlabel=L"\alpha (\epsilon - \bar\epsilon_{\log})",
        ylabel=L"K \geq 2"*" Proportions",
        legend=:bottomright,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm
    )

    # Initialize a new plot where x-axis is alpha(L)*(epsilon - epsilon_log_infinite)
    energy_density_transformed_infinite_plot = scatter(
        xlabel=L"\alpha (L) (\epsilon - \bar\epsilon_{\log}^{\infty})",
        ylabel=L"K \geq 2"*" Proportions",
        legend=:left,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm
    )

    for (i, L) in enumerate(L_values)
        if L == 3
            continue
        end

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




        # If \gamma = fit_K2.param[3] 
        # If \epsilon_log = fit_K2.param[2]
        # If \alpha = fit_K2.param[1]
        alpha = fit_K2.param[1]
        epsilon_log = fit_K2.param[2]
        gamma = fit_K2.param[3]



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
        
        scatter!(collapsed_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            collapsed_K2_values,
            label="L = $L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )

        # Also plot black dashed y=x line 
        plot!(
            collapsed_plot,
            LinRange(-1.0, -0.05, 100),
            LinRange(-1.0, -0.05, 100),
            label="L = $L",
            color=:black,
            linestyle=:dash,
        )



        # TRANSFORMED ENERGY DENSITY PLOT
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


        # INFINITE TRANSFORMED ENERGY DNESITY PLOT 
        alpha_function(L) = 7.9268460128440905 + 2.95038695544134*L
        epsilon_log_infinite = -0.24151854507698275

        transformed_energy_density_infinite = alpha_function(L) * (energy_density_K2 .- epsilon_log_infinite)

        # Plot the transformed energy density against K>=2 proportions
        scatter!(
            energy_density_transformed_infinite_plot,
            transformed_energy_density_infinite[abs.(transformed_energy_density_infinite) .< 1e3],
            k_saddle_proportions[abs.(transformed_energy_density_infinite) .< 1e3, 4],  # Proportion of K >= 2
            label="L = $L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )


   end

    display(collapsed_plot)
    savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_inverted$(addon).pdf")
    savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_inverted$(addon).png")


    display(energy_density_transformed_plot)
    savefig(energy_density_transformed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_transformed_energy_density$(addon).pdf")
    savefig(energy_density_transformed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_transformed_energy_density$(addon).png")

    display(energy_density_transformed_infinite_plot)
    savefig(energy_density_transformed_infinite_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_transformed_energy_density_infinite$(addon).pdf")
    savefig(energy_density_transformed_infinite_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_generalised_logistic_function_transformed_energy_density_infinite$(addon).png")
    
    ##### #####










    ##### COMBINED GENERALISED LOGISTIC FUNCTION PARAMETERS PLOT #####

    # Create a plot for generalised logistic function parameters
    logistic_params_plot = plot(
        layout=(1,3),
        size=(1200, 400),
        margin=5mm,
    )

    # Plot alpha parameters for both K=0 and K≥2 vs L
    plot!(
        logistic_params_plot[1],
        [L for L in L_values],
        [fit_data_dict[L][1].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha",
        title="Alpha Parameter",
        label="K=0",
        marker=:square,
        color=:black,
        legend=:topleft
    )

    # Plot alpha parameter for K≥2 vs L
    plot!(
        logistic_params_plot[1],
        [L for L in L_values],
        [fit_data_dict[L][2].param[1] for L in L_values],
        label="K≥2",
        marker=:diamond,
        color=:blue
    )

    # Plot epsilon_log parameter for K=0 vs 1/L
    plot!(
        logistic_params_plot[2],
        [1/L for L in L_values],
        [fit_data_dict[L][1].param[2] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\bar\epsilon_{\rm log}",
        title="Epsilon Log Parameter",
        label="K=0",
        marker=:square,
        color=:black,
        legend=:bottomleft,
        ylims=(-0.275,-0.1),
        xlims=(0.0,0.22)
    )

    # Plot epsilon_log parameter for K≥2 vs 1/L
    plot!(
        logistic_params_plot[2],
        [1/L for L in L_values],
        [fit_data_dict[L][2].param[2] for L in L_values],
        label="K≥2",
        marker=:diamond,
        color=:blue
    )

    # Plot gamma parameter for K=0 vs 1/L
    plot!(
        logistic_params_plot[3],
        [1/L for L in L_values],
        [fit_data_dict[L][1].param[3] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\gamma",
        title="Gamma Parameter",
        label="K=0",
        marker=:square,
        color=:black,
        legend=:topright,
        xlims=(0.0,0.21)
    )

    # Plot gamma parameter for K≥2 vs 1/L
    plot!(
        logistic_params_plot[3],
        [1/L for L in L_values],
        [fit_data_dict[L][2].param[3] for L in L_values],
        label="K≥2",
        marker=:diamond,
        color=:blue
    )

    # Fit lines to the parameters
    L_vals = [L for L in L_values]
    inverse_L_vals = [1/L for L in L_values]

    # Alpha parameter fits
    alpha_K0_values = [fit_data_dict[L][1].param[1] for L in L_values]
    alpha_K2_values = [fit_data_dict[L][2].param[1] for L in L_values]
    fit_alpha_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_values, [0.0, 0.0])
    fit_alpha_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_values, [0.0, 0.0])

    # Epsilon_log parameter fits
    epsilon_log_K0_values = [fit_data_dict[L][1].param[2] for L in L_values]
    epsilon_log_K2_values = [fit_data_dict[L][2].param[2] for L in L_values]
    fit_epsilon_log_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, epsilon_log_K0_values, [0.0, 0.0])
    fit_epsilon_log_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, epsilon_log_K2_values, [0.0, 0.0])

    # Gamma parameter fits
    gamma_K0_values = [fit_data_dict[L][1].param[3] for L in L_values]
    gamma_K2_values = [fit_data_dict[L][2].param[3] for L in L_values]
    fit_gamma_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K0_values, [0.0, 0.0])
    fit_gamma_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K2_values, [0.0, 0.0])

    # Plot alpha parameter fits
    x_fit_L = LinRange(min(L_values...), max(L_values...), 100)
    plot!(
        logistic_params_plot[1],
        x_fit_L,
        fit_alpha_K0.param[1] .+ fit_alpha_K0.param[2].*x_fit_L,
        label="K=0: $(round(fit_alpha_K0.param[1], digits=2)) + $(round(fit_alpha_K0.param[2], digits=2))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        logistic_params_plot[1],
        x_fit_L,
        fit_alpha_K2.param[1] .+ fit_alpha_K2.param[2].*x_fit_L,
        label="K≥2: $(round(fit_alpha_K2.param[1], digits=2)) + $(round(fit_alpha_K2.param[2], digits=2))·L",
        color=:blue,
        linestyle=:dot
    )

    # Plot epsilon_log and gamma parameter fits
    x_fit_1_L = LinRange(0, 0.22, 100)

    plot!(
        logistic_params_plot[2],
        x_fit_1_L,
        fit_epsilon_log_K0.param[1] .+ fit_epsilon_log_K0.param[2].*x_fit_1_L,
        label="K=0: $(round(fit_epsilon_log_K0.param[1], digits=3)) + $(round(fit_epsilon_log_K0.param[2], digits=3))·(1/L)",
        color=:black,
        linestyle=:dash
    )

    plot!(
        logistic_params_plot[2],
        x_fit_1_L,
        fit_epsilon_log_K2.param[1] .+ fit_epsilon_log_K2.param[2].*x_fit_1_L,
        label="K≥2: $(round(fit_epsilon_log_K2.param[1], digits=3)) + $(round(fit_epsilon_log_K2.param[2], digits=3))·(1/L)",
        color=:blue,
        linestyle=:dot
    )

    plot!(
        logistic_params_plot[3],
        x_fit_1_L,
        fit_gamma_K0.param[1] .+ fit_gamma_K0.param[2].*x_fit_1_L,
        label="K=0: $(round(fit_gamma_K0.param[1], digits=2)) + $(round(fit_gamma_K0.param[2], digits=2))·(1/L)",
        color=:black,
        linestyle=:dash
    )

    plot!(
        logistic_params_plot[3],
        x_fit_1_L,
        fit_gamma_K2.param[1] .+ fit_gamma_K2.param[2].*x_fit_1_L,
        label="K≥2: $(round(fit_gamma_K2.param[1], digits=2)) + $(round(fit_gamma_K2.param[2], digits=2))·(1/L)",
        color=:blue,
        linestyle=:dot
    )

    # Print fit results
    println("K=0 alpha parameter fit: $(round(fit_alpha_K0.param[1], digits=2)) + $(round(fit_alpha_K0.param[2], digits=2))·L")
    println("K≥2 alpha parameter fit: $(round(fit_alpha_K2.param[1], digits=2)) + $(round(fit_alpha_K2.param[2], digits=2))·L")
    println("K=0 epsilon_log parameter fit: $(round(fit_epsilon_log_K0.param[1], digits=3)) + $(round(fit_epsilon_log_K0.param[2], digits=3))·(1/L)")
    println("K≥2 epsilon_log parameter fit: $(round(fit_epsilon_log_K2.param[1], digits=3)) + $(round(fit_epsilon_log_K2.param[2], digits=3))·(1/L)")
    println("K=0 gamma parameter fit: $(round(fit_gamma_K0.param[1], digits=2)) + $(round(fit_gamma_K0.param[2], digits=2))·(1/L)")
    println("K≥2 gamma parameter fit: $(round(fit_gamma_K2.param[1], digits=2)) + $(round(fit_gamma_K2.param[2], digits=2))·(1/L)")

    # Print K>=2 alpha fit parameters (as in original code)
    println("K>=2 alpha fit parameters, p[1] + p[2]L: ", fit_alpha_K2.param)

    # Save and display plots
    savefig(logistic_params_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_generalised_logistic_fit$(addon).pdf")
    savefig(logistic_params_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_generalised_logistic_fit$(addon).png")
    display(logistic_params_plot)
    
    
 



        ##### INSETS ON COMBINED PLOT #####

    # Create the insets first, then add points to them
    # Alpha parameter inset
    scatter!(combined_plot, [], [], 
            inset=bbox(0.25,0.225,0.2,0.2), subplot=2, 
            ylabel=L"\alpha", xlabel=L"L", legend=false, 
            # ylims=(18,57),
            ylims=(10,57), 
            # xlims=(4.5,15.5), 
            xlims=(2.5,15.5), 
            xticks=[3,5,7,9,11,13,15])
    
    # Gamma parameter inset
    scatter!(combined_plot, [], [], 
            inset=bbox(0.75,0.515,0.2,0.2), subplot=3, 
            ylabel=L"\gamma", xlabel=L"1/L", label="", 
            # ylims=(5.0,9.0),
            ylims=(3.0,8.0),
            # xlims=(0.0,0.21),
            xlims=(0.0,0.35), 
            xticks=[0,0.1,0.2,0.3])
    
    # Epsilon_log parameter inset
    scatter!(combined_plot, [], [], 
            inset=bbox(0.255,0.515,0.2,0.2), subplot=4, 
            ylabel=L"\bar\epsilon_{\rm log}", xlabel=L"1/L", label="",
            # ylims=(-0.245,-0.20),
            ylims=(-0.245,-0.165),
            # xlims=(0.0,0.22),
            xlims=(0.0,0.35), 
            xticks=[0,0.1,0.2,0.3])
    
    # Now add points to each inset with the appropriate colors
    # Alpha parameter points
    for (i, L) in enumerate(L_values)
        scatter!(combined_plot, [L], [fit_data_dict[L][2].param[1]], 
                label="", color=colors[i], marker=:diamond, markersize=4, subplot=2)
    end
    
    # Add the linear fit for alpha
    plot!(combined_plot, L_values, fit_alpha_K2.param[1] .+ fit_alpha_K2.param[2].*L_values, 
          label="", color=:black, linestyle=:dash, subplot=2)
    
    # Gamma parameter points
    for (i, L) in enumerate(L_values)
        scatter!(combined_plot, [1/L], [fit_data_dict[L][2].param[3]], 
                label="", color=colors[i], marker=:diamond, markersize=4, subplot=3)
    end
    
    # Epsilon_log parameter points
    for (i, L) in enumerate(L_values)
        scatter!(combined_plot, [1/L], [fit_data_dict[L][2].param[2]], 
                label="", color=colors[i], marker=:diamond, markersize=4, subplot=4)
    end
    
    # Add linear fit for epsilon_log
    epsilon_log_K2_values = [fit_data_dict[L][2].param[2] for L in L_values]
    fit_beta_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, 1 ./ L_values, epsilon_log_K2_values, [0.0, 0.0])
    c, m  = fit_beta_K2.param
    confidence_intervals = confidence_interval(fit_beta_K2, 0.1)
    c_ci = confidence_intervals[2]
    c_error = c - c_ci[1]
    # fitting_1_L_values = LinRange(0, 0.2, 100)
    fitting_1_L_values = LinRange(0, 0.35, 100)
    plot!(combined_plot, fitting_1_L_values, c .+ m.*fitting_1_L_values, label="", color=:black, linestyle=:dash, subplot=4)
    # scatter!(combined_plot, [0.0], [c], yerror=(c_error), subplot=3, color=:blue, markerstrokecolor=:blue, markersize=1, label="", errorbar_color=:blue)
    println("K=2 \epsilon_{log} fit, p[1] + p[2]L: ", fit_beta_K2.param)
    println("L values: ", L_values)
    println("\bar \epsilon_{log} values: ", epsilon_log_K2_values)


    ##### #####


    ##### COMBINED PLOT LEGEND #####

    # Build combined plot legend
    scatter!(combined_plot, [], [], label="K=0", color=:white, marker=:circle, markersize=4)
    scatter!(combined_plot, [], [], label="K≥2", color=:white, marker=:diamond, markersize=4)
    scatter!(combined_plot, [], [], label="L=15", color=colors[1], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=13", color=colors[2], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=11", color=colors[3], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=9", color=colors[4], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=7", color=colors[5], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=5", color=colors[6], marker=:square, markersize=4)
    scatter!(combined_plot, [], [], label="L=3", color=colors[7], marker=:square, markersize=4)


    # Save and display the combined graph
    save_name = "combined_saddle_index_proportions_multiple_L"
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)$(addon).pdf")
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)$(addon).png")
    display(combined_plot)


    ##### #####












    


























    ##### LOG(-EPSILON_LOG) VS LOG(1/L) PLOT #####

    epsilon_log_plot_log = scatter(
    xlabel=L"\log(1/L)",
    ylabel=L"\log(-\bar\epsilon_{\rm log})",
    legend=:bottomright,
    xguidefontsize=12,
    yguidefontsize=12,
    margin=5mm,
    )

    # Plot log(-epsilon_log) values for K=0
    for (i, L) in enumerate(L_values)
        scatter!(epsilon_log_plot_log, [log(1/L)], [log(-fit_data_dict[L][1].param[2])], 
                label=(i==1 ? "K=0" : ""), color=colors[i], marker=:square, markersize=4)
    end

    # Plot log(-epsilon_log) values for K≥2
    for (i, L) in enumerate(L_values)
        scatter!(epsilon_log_plot_log, [log(1/L)], [log(-fit_data_dict[L][2].param[2])], 
                label=(i==1 ? "K≥2" : ""), color=colors[i], marker=:diamond, markersize=4)
    end

    # Add linear fits for log(-epsilon_log)
    epsilon_log_K0_values_log = [log(-fit_data_dict[L][1].param[2]) for L in L_values]
    epsilon_log_K2_values_log = [log(-fit_data_dict[L][2].param[2]) for L in L_values]

    fit_epsilon_log_K0_log = curve_fit((x, p) -> p[1] .+ p[2].*x, log.(1 ./ L_values), epsilon_log_K0_values_log, [0.0, 0.0])
    fit_epsilon_log_K2_log = curve_fit((x, p) -> p[1] .+ p[2].*x, log.(1 ./ L_values), epsilon_log_K2_values_log, [0.0, 0.0])

    fitting_log_1_L_values = LinRange(-3.0, 0, 100)  # Adjusted range for log(1/L)
    plot!(epsilon_log_plot_log, fitting_log_1_L_values, fit_epsilon_log_K0_log.param[1] .+ fit_epsilon_log_K0_log.param[2].*fitting_log_1_L_values, 
        label="K=0: $(round(fit_epsilon_log_K0_log.param[1], digits=3)) + $(round(fit_epsilon_log_K0_log.param[2], digits=3))·log(1/L)", 
        color=:black, linestyle=:dash)
    plot!(epsilon_log_plot_log, fitting_log_1_L_values, fit_epsilon_log_K2_log.param[1] .+ fit_epsilon_log_K2_log.param[2].*fitting_log_1_L_values, 
        label="K≥2: $(round(fit_epsilon_log_K2_log.param[1], digits=3)) + $(round(fit_epsilon_log_K2_log.param[2], digits=3))·log(1/L)", 
        color=:black, linestyle=:dot)

    display(epsilon_log_plot_log)
    savefig(epsilon_log_plot_log, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_generalised_logistic_fit_epsilon_log_log_vs_log_1_L_fixed_epsilon_log_and_gamma$(addon).pdf")
    savefig(epsilon_log_plot_log, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_generalised_logistic_fit_epsilon_log_log_vs_log_1_L_fixed_epsilon_log_and_gamma$(addon).png")


    ##### #####


























    ##### K0 FLOW PLOT #####


    # Choose 10 values between energy density -0.5 and -0.1, and use the fits to plot their value as 1/L changes
    energy_density_values = vcat(LinRange(-0.24-0.2, -0.245, 15), LinRange(-0.25995, -0.23, 3), 0.15)
    L_values = [15,13,11,9,7,5,3]

    # Find the value in energy_density closest to -0.24
    closest_value = energy_density_values[argmin(abs.(energy_density_values .+ 0.24))]


    ## K0 PLOT
    K0_plot = scatter(
        xlabel=L"1/L",
        ylabel=L"K=0"*" Minima Proportions",
            legend=:topleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
    )



    plot!(K0_plot, [], [], color=alex_blue, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon < -0.24")
    for (i, energy_density) in pairs(energy_density_values)
        if energy_density == closest_value
            linestyle = :dash
            color = :black
            label = L"$\epsilon = -0.24$"
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
    plot!(K0_plot, [], [], color=alex_red, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon > -0.24")

    display(K0_plot)
    savefig(K0_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/flow_K0_saddle_index_proportions$(addon).pdf")
    savefig(K0_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/flow_K0_saddle_index_proportions$(addon).png")


    ##### #####








    ##### K2 FLOW PLOT #####

    K2_plot = scatter(
        xlabel=L"1/L",
        ylabel=L"K\geq2"*" Saddle Proportions",
        legend=:topleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
    )

    plot!(K2_plot, [], [], color=alex_red, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon > -0.24")
    for (i, energy_density) in pairs(energy_density_values)
        if energy_density == closest_value
            linestyle = :dash
            color = :black
            label=L"$\epsilon = -0.24$"
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
    plot!(K2_plot, [], [], color=alex_blue, linewidth=3, legend=:right, linestyle=:solid, label=L"\epsilon < -0.24")

    display(K2_plot)
    savefig(K2_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/flow_K2_saddle_index_proportions$(addon).pdf")
    savefig(K2_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/flow_K2_saddle_index_proportions$(addon).png")

    

    ##### #####






















    ##### LOGISTIC GRADIENTS GRAPH #####

    # Calculate the gradient of the logistic function for both K0 and K>=2
    # K0_functional_form(x, p) = (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^(1/p[3])
    # K2_functional_form(x, p) = 1 .- K0_functional_form(x, p)

    # Gradient of K0_functional_form(x, p) = -((p[1]*exp.(p[1] .* (x .- p[2])))/p[3]) .* (1 ./ (1 .+ exp.(p[1] .* (x .- p[2])))).^((1/p[3]) +1)
    # Gradient of K2_functional_form(x, p) = -Gradient of K0_functional_form(x, p)

    # Define analytical gradient functions
    function K0_gradient(x, p)
        # p = [alpha, epsilon_log, gamma]
        alpha, epsilon_log, gamma = p
        # K0 gradient = -((alpha*exp(alpha*(x-epsilon_log)))/gamma) * (1/(1+exp(alpha*(x-epsilon_log))))^(gamma + 1)
        return -((alpha .* exp.(alpha .* (x .- epsilon_log))) .* gamma) .* 
               (1 ./ (1 .+ exp.(alpha .* (x .- epsilon_log)))).^(gamma + 1)
    end

    function K2_gradient(x, p)
        # Since K2 = 1 - K0, the gradient of K2 is the negative of the gradient of K0
        return K0_gradient(x, p)  # We're taking absolute value anyway
    end

    # Create plots for K=0 and K>=2 gradients from analytical form
    K0_gradients_plot = plot(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=0"*" Proportions Gradients",
        legend=:topright,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05),
        ylims=(0, 20)
    )
    
    # Create plot for K>=2 gradients
    K2_gradients_plot = plot(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K\geq 2"*" Proportions Gradients",
        legend=:topright,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05),
        ylims=(0, 20)  # Significantly reduced from 60 to 1 to better show data
    )
    
    L_values = [15, 13, 11, 9, 7, 5]

    # Use the fit parameters to calculate and plot gradients
    gradient_energy_range = LinRange(-0.7, -0.05, 500)  # More points for smoother curves
    
    for (i, L) in enumerate(L_values)
        if !haskey(fit_data_dict, L)
            continue
        end
        
        fit_K0, fit_K2 = fit_data_dict[L]
        
        println("Processing L = $L")
        println("fit_K0.param: ", fit_K0.param)
        println("fit_K2.param: ", fit_K2.param)
        
        # Calculate gradients using analytical expressions
        K0_gradient_values = abs.(K0_gradient(gradient_energy_range, fit_K0.param))
        K2_gradient_values = abs.(K2_gradient(gradient_energy_range, fit_K2.param))
        
        # Plot K=0 gradients
        plot!(
            K0_gradients_plot,
            gradient_energy_range,
            K0_gradient_values,
            label="L = $L",
            color=colors[i],
            linewidth=2
        )
        
        # Find and mark the maximum K=0 gradient point
        max_K0_gradient_idx = argmax(K0_gradient_values)
        max_K0_gradient_energy = gradient_energy_range[max_K0_gradient_idx]
        max_K0_gradient_value = K0_gradient_values[max_K0_gradient_idx]
        
        scatter!(
            K0_gradients_plot,
            [max_K0_gradient_energy],
            [max_K0_gradient_value],
            color=colors[i],
            markersize=6,
            label=""
        )
        
        # Plot K>=2 gradients
        plot!(
            K2_gradients_plot,
            gradient_energy_range,
            K2_gradient_values,
            label="L = $L",
            color=colors[i],
            linewidth=2
        )
        
        # Find and mark the maximum K>=2 gradient point
        max_K2_gradient_idx = argmax(K2_gradient_values)
        max_K2_gradient_energy = gradient_energy_range[max_K2_gradient_idx]
        max_K2_gradient_value = K2_gradient_values[max_K2_gradient_idx]
        
        scatter!(
            K2_gradients_plot,
            [max_K2_gradient_energy],
            [max_K2_gradient_value],
            color=colors[i],
            markersize=6,
            label=""
        )
        
        # Print the maximum gradient points
        println("L=$L: Maximum K=0 gradient at energy density ϵ = $(round(max_K0_gradient_energy, digits=3)) with value $(round(max_K0_gradient_value, digits=3))")
        println("L=$L: Maximum K>=2 gradient at energy density ϵ = $(round(max_K2_gradient_energy, digits=3)) with value $(round(max_K2_gradient_value, digits=3))")
    end

    # Add vertical line at ϵ = -0.24 to both plots
    vline!(K0_gradients_plot, [-0.24], linecolor=:black, linestyle=:dash, linewidth=1, label=L"\epsilon = -0.24")
    vline!(K2_gradients_plot, [-0.24], linecolor=:black, linestyle=:dash, linewidth=1, label=L"\epsilon = -0.24")
    
    # Save and display the K=0 gradients plot
    savefig(K0_gradients_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/gradients_K0_generalised_logistic_function_by_L$(addon).pdf")
    savefig(K0_gradients_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/gradients_K0_generalised_logistic_function_by_L$(addon).png")
    display(K0_gradients_plot)
    
    # Save and display the K>=2 gradients plot
    savefig(K2_gradients_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/gradients_K2_generalised_logistic_function_by_L$(addon).pdf")
    savefig(K2_gradients_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/gradients_K2_generalised_logistic_function_by_L$(addon).png")
    display(K2_gradients_plot)


    ##### #####


























































    ##### Erf(x) FIT TO K0 and K2 #####
    # Use K2_functional_form(x, p) = 0.5 + 0.5*erf(p[1]*(x-p[2]))
    # Use K0_functional_form(x, p) = 1 - K2_functional_form(x, p)

    # Define erf-based functional forms
    erf_K2_functional_form(x, p) = 0.5 .+ 0.5 .* erf.(p[1] .* (x .- p[2]))
    erf_K0_functional_form(x, p) = 0.5 .- 0.5 .* erf.(p[1] .* (x .- p[2]))  # Note: This is NOT 1 - K2 form, it's a separate fit

    # Initial parameters [scale, midpoint]
    p0_erf = [10.0, -0.24]
    
    # Create a combined plot for erf fits
    erf_combined_plot = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
        legend=:left,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05)
    )
    
    # Dictionary to store erf fit parameters
    erf_fit_data_dict = Dict()
    
    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Fit erf functions to K≥2 data
        valid_rows_K2 = .!isnan.(k_saddle_proportions[:, 4]) .& .!isinf.(k_saddle_proportions[:, 4]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K2 = E0_bin_values[valid_rows_K2]
        k_saddle_proportions_filtered_K2 = k_saddle_proportions[valid_rows_K2, 4]
        
        # Transform E0_bin_values to energy density
        energy_density_K2 = -1.0 .+ (E0_bin_values_filtered_K2./-solved_configuration_energy(RubiksCube(L)))
        
        # Fit K≥2 data with erf function
        erf_fit_K2 = curve_fit(erf_K2_functional_form, energy_density_K2, k_saddle_proportions_filtered_K2, p0_erf, show_trace=false, maxIter=200)
        
        # Fit erf functions to K0 data
        valid_rows_K0 = .!isnan.(k_saddle_proportions[:, 1]) .& .!isinf.(k_saddle_proportions[:, 1]) .& .!isnan.(E0_bin_values) .& .!isinf.(E0_bin_values)
        E0_bin_values_filtered_K0 = E0_bin_values[valid_rows_K0]
        k_saddle_proportions_filtered_K0 = k_saddle_proportions[valid_rows_K0, 1]
        
        # Transform E0_bin_values to energy density
        energy_density_K0 = -1.0 .+ (E0_bin_values_filtered_K0./-solved_configuration_energy(RubiksCube(L)))
        
        # Fit K0 data with erf function
        erf_fit_K0 = curve_fit(erf_K0_functional_form, energy_density_K0, k_saddle_proportions_filtered_K0, p0_erf, show_trace=false, maxIter=200)
        
        # Store both fit parameters
        erf_fit_data_dict[L] = (erf_fit_K0, erf_fit_K2)
        
        # Plot the data points
        scatter!(erf_combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 1],
            label="L = $L",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
        
        scatter!(erf_combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            label="L = $L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )

        # Plot fitted curves
        fit_values = LinRange(-1.0, -0.05, 100)
        plot!(
                erf_combined_plot,
                fit_values,
                erf_K0_functional_form(fit_values, erf_fit_K0.param),
                label=i==1 ? "K=0 fit" : "",
                color=colors[i],
                linestyle=:dash,
            )
            
        plot!(
                erf_combined_plot,
                fit_values,
                erf_K2_functional_form(fit_values, erf_fit_K2.param),
                label=i==1 ? "K≥2 fit" : "",
                color=colors[i],
                linestyle=:dash,
            )



        
        # Find the midpoints (where K0 = 0.5 and K2 = 0.5)
        midpoint_K0 = erf_fit_K0.param[2]
        midpoint_K2 = erf_fit_K2.param[2]
        println("L=$L: K0 erf midpoint at energy density ϵ = $(round(midpoint_K0, digits=4))")
        println("L=$L: K2 erf midpoint at energy density ϵ = $(round(midpoint_K2, digits=4))")
        
        # Find the scaling factors
        scale_K0 = erf_fit_K0.param[1]
        scale_K2 = erf_fit_K2.param[1]
        println("L=$L: K0 erf scale factor = $(round(scale_K0, digits=4))")
        println("L=$L: K2 erf scale factor = $(round(scale_K2, digits=4))")
        
        # Create individual plot for this L value showing both K0 and K2 data and erf fits
        individual_erf_plot = scatter(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
            ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
            title="L = $L",
            legend=:topleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05),
            ylims=(0, 1.05)
        )
        
        # Plot the K0 data points
        scatter!(individual_erf_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 1],
            label="K=0 data",
            color=colors[i],
            marker=:circle,
            markersize=4,
        )
        
        # Plot the K2 data points
        scatter!(individual_erf_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            label="K≥2 data",
            color=colors[i],
        marker=:diamond,
            markersize=4,
    )

        # Plot the K0 erf fit with black dashed line
    plot!(
            individual_erf_plot,
            fit_values,
            erf_K0_functional_form(fit_values, erf_fit_K0.param),
            label="K=0 erf fit: 0.5 - 0.5·erf($(round(scale_K0, digits=2))·(ϵ-$(round(midpoint_K0, digits=4))))",
        color=:black,
            linestyle=:dash,
            linewidth=2
    )

        # Plot the K2 erf fit with blue dashed line
    plot!(
            individual_erf_plot,
            fit_values,
            erf_K2_functional_form(fit_values, erf_fit_K2.param),
            label="K≥2 erf fit: 0.5 + 0.5·erf($(round(scale_K2, digits=2))·(ϵ-$(round(midpoint_K2, digits=4))))",
            color=:blue,
            linestyle=:dash,
            linewidth=2
        )
        
        # Add vertical lines at the midpoints
        vline!(individual_erf_plot, [midpoint_K0], linecolor=:black, linestyle=:dot, linewidth=1, 
               label="K=0 midpoint: ϵ = $(round(midpoint_K0, digits=4))")
        vline!(individual_erf_plot, [midpoint_K2], linecolor=:blue, linestyle=:dot, linewidth=1, 
               label="K≥2 midpoint: ϵ = $(round(midpoint_K2, digits=4))")
        
        # Add horizontal line at 0.5
        hline!(individual_erf_plot, [0.5], linecolor=:gray, linestyle=:dot, linewidth=1, label="")
        
        # Save the individual plot
        savefig(individual_erf_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_erf_fits$(addon).pdf")
        savefig(individual_erf_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_erf_fits$(addon).png")
        display(individual_erf_plot)
    end

    # Create a plot for scale and midpoint parameters versus L
    erf_params_plot = plot(
        layout=(2,2),
        size=(900, 800),
        margin=5mm,
    )

    # Plot K0 midpoint vs 1/L
    plot!(
        erf_params_plot[1],
        [1/L for L in L_values],
        [erf_fit_data_dict[L][1].param[2] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\epsilon_{\text{mid}}^{K=0}",
        title="K=0 Midpoint Parameter",
        label="Data",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )

    # Plot K2 midpoint vs 1/L
    plot!(
        erf_params_plot[2],
        [1/L for L in L_values],
        [erf_fit_data_dict[L][2].param[2] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\epsilon_{\text{mid}}^{K\geq2}",
        title="K≥2 Midpoint Parameter",
        label="Data",
        marker=:circle,
        color=:blue,
        legend=:bottomright
    )

    # Plot K0 scale factor vs L
    plot!(
        erf_params_plot[3],
        [L for L in L_values],
        [erf_fit_data_dict[L][1].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha^{K=0}",
        title="K=0 Scale Parameter",
        label="Data",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )

    # Plot K2 scale factor vs L
    plot!(
        erf_params_plot[4],
        [L for L in L_values],
        [erf_fit_data_dict[L][2].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha^{K\geq2}",
        title="K≥2 Scale Parameter",
        label="Data",
        marker=:circle,
        color=:blue,
        legend=:bottomright
    )

    # Fit lines to the parameters
    L_vals = [L for L in L_values]
    inverse_L_vals = [1/L for L in L_values]

    # K0 scale factor fit
    scale_K0_vals = [erf_fit_data_dict[L][1].param[1] for L in L_values]
    fit_scale_K0 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, scale_K0_vals, [0.0, 0.0])

    # K2 scale factor fit
    scale_K2_vals = [erf_fit_data_dict[L][2].param[1] for L in L_values]
    fit_scale_K2 = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, scale_K2_vals, [0.0, 0.0])

    # Plot scale factor fits
    x_fit_L = LinRange(min(L_values...), max(L_values...), 100)
    plot!(
        erf_params_plot[3],
        x_fit_L,
        fit_scale_K0.param[1] .+ fit_scale_K0.param[2].*x_fit_L,
        label="$(round(fit_scale_K0.param[1], digits=4)) + $(round(fit_scale_K0.param[2], digits=4))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        erf_params_plot[4],
        x_fit_L,
        fit_scale_K2.param[1] .+ fit_scale_K2.param[2].*x_fit_L,
        label="$(round(fit_scale_K2.param[1], digits=4)) + $(round(fit_scale_K2.param[2], digits=4))·L",
        color=:blue,
        linestyle=:dash
    )


    # Print fit results
    println("K=0 scale parameter fit: $(round(fit_scale_K0.param[1], digits=4)) + $(round(fit_scale_K0.param[2], digits=4))·L")
    println("K≥2 scale parameter fit: $(round(fit_scale_K2.param[1], digits=4)) + $(round(fit_scale_K2.param[2], digits=4))·L")



    # Save and display plots
    savefig(erf_combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_erf_fit$(addon).pdf")
    savefig(erf_combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_saddle_index_proportions_erf_fit$(addon).png")
    display(erf_combined_plot)

    savefig(erf_params_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_erf_fit$(addon).pdf")
    savefig(erf_params_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_erf_fit$(addon).png")
    display(erf_params_plot)
    
 
    ##### COLLAPSED PLOT OF K2 PROPORTIONS USING ERF SCALE FACTOR AND MIDPOINT #####

    # Create collapsed plot where x-axis is scale_factor*(epsilon - midpoint)
    erf_collapsed_plot = scatter(
        xlabel=L"\alpha \cdot (\epsilon - \epsilon_{\rm{mid}})",
        ylabel=L"K \geq 2"*" Proportions",
        legend=:bottomright,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-15,)
    )

    for (i, L) in enumerate(L_values)
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Get the erf fit parameters for K≥2
        erf_fit_K0, erf_fit_K2 = erf_fit_data_dict[L]
        scale_factor = erf_fit_K2.param[1]
        midpoint = erf_fit_K2.param[2]
        
        # Transform E0_bin_values to energy density
        energy_density = -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L)))
        
        # Calculate the collapsed x-axis values: scale_factor * (epsilon - midpoint)
        collapsed_x_values = scale_factor .* (energy_density .- midpoint)

        valid_indices = abs.(collapsed_x_values) .< 15
        
        # Plot the collapsed data
        scatter!(erf_collapsed_plot,
            collapsed_x_values[valid_indices],
            k_saddle_proportions[:, 4],  # K≥2 proportions
            label="L = $L",
            color=colors[i],
            marker=:diamond,
            markersize=4,
            alpha=0.7
        )
    end

    # Save and display the collapsed plot
    savefig(erf_collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_erf_K2_proportions$(addon).pdf")
    savefig(erf_collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_erf_K2_proportions$(addon).png")
    display(erf_collapsed_plot)

    println("Erf collapsed plot created - all K≥2 data should collapse onto the universal erf curve")

    ##### #####



























    


























    ##### F=0.5 GRADIENT ESTIMATIONS OF EPSILON_LOG FOR BOTH K0 AND K2 #####

    # Function to estimate epsilon_log using gradients
    function estimate_epsilon_log_from_gradient(func, fit_params, fit_name, L, color_idx, addon)
        # Create a plot for this fit
        gradient_estimate_plot = plot(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
            ylabel=L"$fit_name"*" Proportions",
            title="L = $L, $fit_name",
            legend=:topleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05),
            ylims=(0, 1.05)
        )
        
        # Plot the function
        energy_range = LinRange(-0.7, -0.05, 500)
        plot!(
            gradient_estimate_plot,
            energy_range,
            func(energy_range, fit_params),
            label="$fit_name fit",
            color=colors[color_idx],
            linewidth=2,
            linestyle=:dash
        )

        # Add the data points
        E0_bin_values, k_saddle_proportions = data_dict[L]
        data_idx = fit_name == "K=0" ? 1 : 4
        scatter!(
            gradient_estimate_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, data_idx],
            label="$fit_name data",
            color=colors[color_idx],
        )
        
        # Find the energy density where function = 0.5
        half_difference(x) = func(x, fit_params) - 0.5
        epsilon_half = find_zero(half_difference, (-0.7, -0.05), Bisection())
        
        # Get parameters for gradient calculation
        alpha, epsilon_log, gamma = fit_params
        
        # Calculate gradient analytically
        # K0: gradient = -(alpha*gamma*exp(alpha*(x-epsilon_log)))/(1+exp(alpha*(x-epsilon_log)))^(gamma+1)
        # K2: gradient = (alpha*gamma*exp(alpha*(x-epsilon_log)))/(1+exp(alpha*(x-epsilon_log)))^(gamma+1)
        # Sign depends on the function (negative for K0, positive for K2)
        sign_factor = fit_name == "K=0" ? -1 : 1
        gradient_at_half = sign_factor * (alpha * gamma * exp(alpha * (epsilon_half - epsilon_log))) / 
                          ((1 + exp(alpha * (epsilon_half - epsilon_log)))^(gamma + 1))
        
        # Calculate the tangent line at function = 0.5
        tangent_line(x) = 0.5 + gradient_at_half * (x - epsilon_half)
        
        # Find where the tangent line intersects function = target value (0 for K0, 1 for K2)
        target_val = fit_name == "K=0" ? 0.0 : 1.0
        epsilon_log_estimate = epsilon_half + (target_val - 0.5) / gradient_at_half
        
        # Plot the tangent line
        tangent_range = LinRange(epsilon_half - 0.2, epsilon_log_estimate + 0.05, 100)
        plot!(
            gradient_estimate_plot,
            tangent_range,
            tangent_line.(tangent_range),
            label="Tangent at $fit_name = 0.5",
            color=:black,
            linestyle=:dash,
            linewidth=2
        )
        
        # Mark the function = 0.5 point
        scatter!(
            gradient_estimate_plot,
            [epsilon_half],
            [0.5],
            label="$fit_name = 0.5 at ϵ = $(round(epsilon_half, digits=4))",
            color=:red,
            markersize=6
        )
        
        # Mark the epsilon_log estimate
        scatter!(
            gradient_estimate_plot,
            [epsilon_log_estimate],
            [target_val],
            label=L"\bar\epsilon_{\rm log}^{\rm est} = $(round(epsilon_log_estimate, digits=4))",
            color=:blue,
            markersize=6
        )
        
        # Add horizontal lines
        hline!(gradient_estimate_plot, [0.5], linecolor=:gray, linestyle=:dot, linewidth=1, label="")
        hline!(gradient_estimate_plot, [target_val], linecolor=:gray, linestyle=:dot, linewidth=1, label="")
        
        # Add the actual epsilon_log from the fit
        vline!(gradient_estimate_plot, [fit_params[2]], linecolor=:green, linestyle=:dash, linewidth=1, 
               label=L"\bar\epsilon_{\rm log} = $(round(fit_params[2], digits=4))")
        
        # Save and display the plot
        file_prefix = fit_name == "K=0" ? "K0" : "K2"
        savefig(gradient_estimate_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(file_prefix)_epsilon_log_gradient_estimate$(addon).pdf")
        savefig(gradient_estimate_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(file_prefix)_epsilon_log_gradient_estimate$(addon).png")
        display(gradient_estimate_plot)
        
        # Print the results
        println("L=$L: $fit_name = 0.5 at ϵ = $(round(epsilon_half, digits=4))")
        println("L=$L: Gradient at $fit_name = 0.5: $(round(gradient_at_half, digits=4))")
        println("L=$L: ϵ_log estimate = $(round(epsilon_log_estimate, digits=4))")
        println("L=$L: Actual ϵ_log from fit = $(round(fit_params[2], digits=4))")
        
        return epsilon_log_estimate
    end
    
    
    # Dictionaries to store epsilon_log estimates
    epsilon_log_K0_estimates = Dict()
    epsilon_log_K2_estimates = Dict()
    
    for (i, L) in enumerate(L_values)
        # Get the fit parameters
        fit_K0, fit_K2 = fit_data_dict[L]
        
        # Estimate epsilon_log for K0
        epsilon_log_K0_estimate = estimate_epsilon_log_from_gradient(
            K0_functional_form, fit_K0.param, "K=0", L, i, addon)
        epsilon_log_K0_estimates[L] = epsilon_log_K0_estimate
        
        # Estimate epsilon_log for K2
        epsilon_log_K2_estimate = estimate_epsilon_log_from_gradient(
            K2_functional_form, fit_K2.param, "K≥2", L, i, addon)
        epsilon_log_K2_estimates[L] = epsilon_log_K2_estimate
        
    end
    
    # Process K0 estimates
    inverse_L_vals = [1/L for L in L_values]
    epsilon_log_K0_estimate_values = [epsilon_log_K0_estimates[L] for L in L_values]
    
    # Process K2 estimates
    epsilon_log_K2_estimate_values = [epsilon_log_K2_estimates[L] for L in L_values]
    





















    

    ##### EPSILON_LOG FIXED LOGISTIC FIT FOR BOTH K0 AND K2 #####

    # Define fixed epsilon log function for K0 (opposite sign compared to K2)
    function K0_fixed_epsilon_log(x, p, epsilon_log_fixed)
        # p[1] = alpha, p[2] = gamma
        return (1 ./ (1 .+ exp.(p[1] .* (x .- epsilon_log_fixed)))).^(p[2])
    end

    # Dictionaries to store the new fits with fixed epsilon_log
    fixed_epsilon_log_K0_fits = Dict()
    fixed_epsilon_log_K2_fits = Dict()
    
    # Function to fit data with fixed epsilon_log and create comparison plots
    function fit_with_fixed_epsilon_log(L, i, fit_params, fit_func, fixed_epsilon_log_func, 
                                        epsilon_log_estimate, saddle_index, fit_name, addon)
        # Get the data for this L value
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Filter valid data
        valid_rows = .!isnan.(k_saddle_proportions[:, saddle_index]) .& 
                    .!isinf.(k_saddle_proportions[:, saddle_index]) .& 
                    .!isnan.(E0_bin_values) .& 
                    .!isinf.(E0_bin_values)
        E0_bin_values_filtered = E0_bin_values[valid_rows]
        k_saddle_proportions_filtered = k_saddle_proportions[valid_rows, saddle_index]
        
        # Transform E0_bin_values to energy density
        energy_density = -1.0 .+ (E0_bin_values_filtered./-solved_configuration_energy(RubiksCube(L)))
        
        # Initial guess for [alpha, gamma]
        p0_fixed = [fit_params[1], fit_params[3]]
        
        # Fit the data with fixed epsilon_log
        fixed_fit = curve_fit((x, p) -> fixed_epsilon_log_func(x, p, epsilon_log_estimate), 
                               energy_density, k_saddle_proportions_filtered, 
                               p0_fixed, show_trace=false, maxIter=200)
        
        # Create a plot comparing the original and fixed fits
        fixed_fit_plot = plot(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
            ylabel="$fit_name Proportions",
            title="L = $L, $fit_name",
            legend=:topleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05),
            ylims=(0, 1.05)
        )
        
        # Plot the data points
        scatter!(fixed_fit_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, saddle_index],
            label="Data",
            color=:gray,
            marker=:circle,
            markersize=4,
            alpha=0.7
        )
        
        # Plot the original fit
        energy_range = LinRange(-0.7, -0.05, 500)
        plot!(
            fixed_fit_plot,
            energy_range,
            fit_func(energy_range, fit_params),
            label="Original fit (ϵ_log=$(round(fit_params[2], digits=4)), γ=$(round(fit_params[3], digits=2)))",
            color=:black,
            linewidth=2
        )
        
        # Plot the fixed epsilon_log fit
        plot!(
            fixed_fit_plot,
            energy_range,
            fixed_epsilon_log_func(energy_range, fixed_fit.param, epsilon_log_estimate),
            label="Fixed ϵ_log fit (ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(fixed_fit.param[2], digits=2)))",
            color=:blue,
            linewidth=2,
            linestyle=:dash
        )
        
        # Add vertical lines for the two epsilon_log values
        vline!(fixed_fit_plot, [fit_params[2]], linecolor=:black, linestyle=:dash, linewidth=1, 
               label="Original ϵ_log")
        vline!(fixed_fit_plot, [epsilon_log_estimate], linecolor=:blue, linestyle=:dash, linewidth=1, 
               label="Fixed ϵ_log")
        
        # Calculate and display the fit quality metrics
        original_rmse = sqrt(sum((fit_func(energy_density, fit_params) .- k_saddle_proportions_filtered).^2) / 
                              length(k_saddle_proportions_filtered))
        fixed_rmse = sqrt(sum((fixed_epsilon_log_func(energy_density, fixed_fit.param, epsilon_log_estimate) .- 
                                k_saddle_proportions_filtered).^2) / length(k_saddle_proportions_filtered))
        
        # Print the fit parameters and quality metrics
        println("L=$L: Original $fit_name fit: α=$(round(fit_params[1], digits=2)), ϵ_log=$(round(fit_params[2], digits=4)), γ=$(round(fit_params[3], digits=2))")
        println("L=$L: Fixed $fit_name fit: α=$(round(fixed_fit.param[1], digits=2)), ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(fixed_fit.param[2], digits=2))")
        println("L=$L: Original RMSE: $(round(original_rmse, digits=5))")
        println("L=$L: Fixed RMSE: $(round(fixed_rmse, digits=5))")
        println("L=$L: RMSE difference: $(round(fixed_rmse - original_rmse, digits=5))")
        
        # Save and display the plot
        savefig(fixed_fit_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(fit_name)_fixed_epsilon_log_fit$(addon).pdf")
        savefig(fixed_fit_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(fit_name)_fixed_epsilon_log_fit$(addon).png")
        display(fixed_fit_plot)
        
        return fixed_fit
    end
    
    # Fit both K0 and K2 data with fixed epsilon_log
    for (i, L) in enumerate(L_values)
        # Get the original fit parameters
        fit_K0, fit_K2 = fit_data_dict[L]
        
        # Get the epsilon_log estimates
        epsilon_log_K0_estimate = epsilon_log_K0_estimates[L]
        epsilon_log_K2_estimate = epsilon_log_K2_estimates[L]
        
        # Fit K0 data with fixed epsilon_log
        fixed_fit_K0 = fit_with_fixed_epsilon_log(
            L, i, fit_K0.param, K0_functional_form, K0_fixed_epsilon_log,
            epsilon_log_K0_estimate, 1, "K=0", addon
        )
        fixed_epsilon_log_K0_fits[L] = fixed_fit_K0
        
        # Fit K2 data with fixed epsilon_log
        fixed_fit_K2 = fit_with_fixed_epsilon_log(
            L, i, fit_K2.param, K2_functional_form, K2_fixed_epsilon_log,
            epsilon_log_K2_estimate, 4, "K≥2", addon
        )
        fixed_epsilon_log_K2_fits[L] = fixed_fit_K2
    end
    
    # Create a plot comparing the alpha and gamma parameters from the original and fixed fits for both K0 and K2
    params_comparison_plot = plot(
        layout=(2,2),
        size=(900, 800),
        margin=5mm,
    )

    # Plot alpha vs L for K0
    plot!(
        params_comparison_plot[1],
        [L for L in L_values],
        [fit_data_dict[L][1].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha \textrm{ for } K=0",
        title="Alpha Parameter (K=0)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )

    plot!(
        params_comparison_plot[1],
        [L for L in L_values],
        [fixed_epsilon_log_K0_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )
    
    # Plot alpha vs L for K2
    plot!(
        params_comparison_plot[2],
        [L for L in L_values],
        [fit_data_dict[L][2].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha \textrm{ for } K \geq 2",
        title="Alpha Parameter (K≥2)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )
    
    plot!(
        params_comparison_plot[2],
        [L for L in L_values],
        [fixed_epsilon_log_K2_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    # Plot gamma vs 1/L for K0
    plot!(
        params_comparison_plot[3],
        [1/L for L in L_values],
        [fit_data_dict[L][1].param[3] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\gamma \textrm{ for } K=0",
        title="Gamma Parameter (K=0)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:topright
    )

    plot!(
        params_comparison_plot[3],
        [1/L for L in L_values],
        [fixed_epsilon_log_K0_fits[L].param[2] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    # Plot gamma vs 1/L for K2
    plot!(
        params_comparison_plot[4],
        [1/L for L in L_values],
        [fit_data_dict[L][2].param[3] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\gamma \textrm{ for } K \geq 2",
        title="Gamma Parameter (K≥2)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:topright
    )

    plot!(
        params_comparison_plot[4],
        [1/L for L in L_values],
        [fixed_epsilon_log_K2_fits[L].param[2] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    # Fit lines to the parameters
    L_vals = [L for L in L_values]
    inverse_L_vals = [1/L for L in L_values]

    # K0 fits
    alpha_K0_orig_vals = [fit_data_dict[L][1].param[1] for L in L_values]
    fit_alpha_K0_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_orig_vals, [0.0, 0.0])
    
    alpha_K0_fixed_vals = [fixed_epsilon_log_K0_fits[L].param[1] for L in L_values]
    fit_alpha_K0_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_fixed_vals, [0.0, 0.0])
    
    gamma_K0_orig_vals = [fit_data_dict[L][1].param[3] for L in L_values]
    fit_gamma_K0_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K0_orig_vals, [0.0, 0.0])
    
    gamma_K0_fixed_vals = [fixed_epsilon_log_K0_fits[L].param[2] for L in L_values]
    fit_gamma_K0_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K0_fixed_vals, [0.0, 0.0])
    
    # K2 fits
    alpha_K2_orig_vals = [fit_data_dict[L][2].param[1] for L in L_values]
    fit_alpha_K2_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_orig_vals, [0.0, 0.0])
    
    alpha_K2_fixed_vals = [fixed_epsilon_log_K2_fits[L].param[1] for L in L_values]
    fit_alpha_K2_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_fixed_vals, [0.0, 0.0])
    
    gamma_K2_orig_vals = [fit_data_dict[L][2].param[3] for L in L_values]
    fit_gamma_K2_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K2_orig_vals, [0.0, 0.0])
    
    gamma_K2_fixed_vals = [fixed_epsilon_log_K2_fits[L].param[2] for L in L_values]
    fit_gamma_K2_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, inverse_L_vals, gamma_K2_fixed_vals, [0.0, 0.0])
    
    # Plot fit lines
    x_fit_L = LinRange(min(L_values...), max(L_values...), 100)
    
    # K0 alpha fit lines
    plot!(
        params_comparison_plot[1],
        x_fit_L,
        fit_alpha_K0_orig.param[1] .+ fit_alpha_K0_orig.param[2].*x_fit_L,
        label="Original: $(round(fit_alpha_K0_orig.param[1], digits=2)) + $(round(fit_alpha_K0_orig.param[2], digits=2))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        params_comparison_plot[1],
        x_fit_L,
        fit_alpha_K0_fixed.param[1] .+ fit_alpha_K0_fixed.param[2].*x_fit_L,
        label="Fixed: $(round(fit_alpha_K0_fixed.param[1], digits=2)) + $(round(fit_alpha_K0_fixed.param[2], digits=2))·L",
        color=:blue,
        linestyle=:dash
    )
    
    # K2 alpha fit lines
    plot!(
        params_comparison_plot[2],
        x_fit_L,
        fit_alpha_K2_orig.param[1] .+ fit_alpha_K2_orig.param[2].*x_fit_L,
        label="Original: $(round(fit_alpha_K2_orig.param[1], digits=2)) + $(round(fit_alpha_K2_orig.param[2], digits=2))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        params_comparison_plot[2],
        x_fit_L,
        fit_alpha_K2_fixed.param[1] .+ fit_alpha_K2_fixed.param[2].*x_fit_L,
        label="Fixed: $(round(fit_alpha_K2_fixed.param[1], digits=2)) + $(round(fit_alpha_K2_fixed.param[2], digits=2))·L",
        color=:blue,
        linestyle=:dash
    )
    
    # Gamma fit lines
    x_fit = LinRange(0, 0.22, 100)
    
    # K0 gamma fit lines
    plot!(
        params_comparison_plot[3],
        x_fit,
        fit_gamma_K0_orig.param[1] .+ fit_gamma_K0_orig.param[2].*x_fit,
        label="Original: $(round(fit_gamma_K0_orig.param[1], digits=2)) + $(round(fit_gamma_K0_orig.param[2], digits=2))·(1/L)",
        color=:black,
        linestyle=:dash
    )

    plot!(
        params_comparison_plot[3],
        x_fit,
        fit_gamma_K0_fixed.param[1] .+ fit_gamma_K0_fixed.param[2].*x_fit,
        label="Fixed: $(round(fit_gamma_K0_fixed.param[1], digits=2)) + $(round(fit_gamma_K0_fixed.param[2], digits=2))·(1/L)",
        color=:blue,
        linestyle=:dash
    )

    # K2 gamma fit lines
    plot!(
        params_comparison_plot[4],
        x_fit,
        fit_gamma_K2_orig.param[1] .+ fit_gamma_K2_orig.param[2].*x_fit,
        label="Original: $(round(fit_gamma_K2_orig.param[1], digits=2)) + $(round(fit_gamma_K2_orig.param[2], digits=2))·(1/L)",
        color=:black,
        linestyle=:dash
    )

    plot!(
        params_comparison_plot[4],
        x_fit,
        fit_gamma_K2_fixed.param[1] .+ fit_gamma_K2_fixed.param[2].*x_fit,
        label="Fixed: $(round(fit_gamma_K2_fixed.param[1], digits=2)) + $(round(fit_gamma_K2_fixed.param[2], digits=2))·(1/L)",
        color=:blue,
        linestyle=:dash
    )

    # Print fit results
    println("K0 original alpha fit: $(round(fit_alpha_K0_orig.param[1], digits=4)) + $(round(fit_alpha_K0_orig.param[2], digits=4))·L")
    println("K0 fixed alpha fit: $(round(fit_alpha_K0_fixed.param[1], digits=4)) + $(round(fit_alpha_K0_fixed.param[2], digits=4))·L")
    println("K0 original gamma fit: $(round(fit_gamma_K0_orig.param[1], digits=4)) + $(round(fit_gamma_K0_orig.param[2], digits=4))·(1/L)")
    println("K0 fixed gamma fit: $(round(fit_gamma_K0_fixed.param[1], digits=4)) + $(round(fit_gamma_K0_fixed.param[2], digits=4))·(1/L)")
    
    println("K2 original alpha fit: $(round(fit_alpha_K2_orig.param[1], digits=4)) + $(round(fit_alpha_K2_orig.param[2], digits=4))·L")
    println("K2 fixed alpha fit: $(round(fit_alpha_K2_fixed.param[1], digits=4)) + $(round(fit_alpha_K2_fixed.param[2], digits=4))·L")
    println("K2 original gamma fit: $(round(fit_gamma_K2_orig.param[1], digits=4)) + $(round(fit_gamma_K2_orig.param[2], digits=4))·(1/L)")
    println("K2 fixed gamma fit: $(round(fit_gamma_K2_fixed.param[1], digits=4)) + $(round(fit_gamma_K2_fixed.param[2], digits=4))·(1/L)")
    
    # Save and display the comparison plot
    savefig(params_comparison_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_fixed_vs_original$(addon).pdf")
    savefig(params_comparison_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_fixed_vs_original$(addon).png")
    display(params_comparison_plot)
    

    ##### #####













    ##### EPSILON_LOG AND GAMMA FIXED LOGISTIC FIT FOR BOTH K0 AND K2 #####

    # Define functions for fixed epsilon_log and gamma fits
    function K0_fixed_epsilon_log_and_gamma(x, p, epsilon_log_fixed, gamma_fixed)
        # p[1] = alpha
        return (1 ./ (1 .+ exp.(p[1] .* (x .- epsilon_log_fixed)))).^(gamma_fixed)
    end

    # Do another fit of both K0 and K2 proportions using the logistic function but with BOTH the epsilon_log parameter 
    # fixed to the estimate from the gradient and the gamma parameter fixed to the mean gamma value from the fixed epsilon_log fits
    # Plot new graphs for each L with the fixed epsilon_log and fixed gamma fit in red, the fixed epsilon_log fit in blue and the original fit in black





    # Dictionaries to store the new fits with fixed epsilon_log and gamma
    fixed_epsilon_log_and_gamma_K0_fits = Dict()
    fixed_epsilon_log_and_gamma_K2_fits = Dict()
    
    # Function to fit with fixed epsilon_log and gamma and create comparison plots
    function fit_with_fixed_epsilon_log_and_gamma(L, i, fit_params, fit_func, fixed_epsilon_log_func,
                                                fixed_epsilon_log_and_gamma_func, epsilon_log_estimate, 
                                                mean_gamma, saddle_index, fit_name, addon)
        # Get the data for this L value
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Filter valid data
        valid_rows = .!isnan.(k_saddle_proportions[:, saddle_index]) .& 
                    .!isinf.(k_saddle_proportions[:, saddle_index]) .& 
                    .!isnan.(E0_bin_values) .& 
                    .!isinf.(E0_bin_values)
        E0_bin_values_filtered = E0_bin_values[valid_rows]
        k_saddle_proportions_filtered = k_saddle_proportions[valid_rows, saddle_index]
        
        # Transform E0_bin_values to energy density
        energy_density = -1.0 .+ (E0_bin_values_filtered./-solved_configuration_energy(RubiksCube(L)))
        
        # Initial guess for [alpha]
        p0_fixed = [fit_params[1]]
        
        # Fit the data with fixed epsilon_log and gamma
        fixed_fit = curve_fit((x, p) -> fixed_epsilon_log_and_gamma_func(x, p, epsilon_log_estimate, mean_gamma), 
                             energy_density, k_saddle_proportions_filtered, 
                             p0_fixed, show_trace=false, maxIter=200)
        
        return fixed_fit
    end
        
    # Apply the fits to both K0 and K2 data
    for (i, L) in enumerate(L_values)
        # Get the original fit parameters
        fit_K0, fit_K2 = fit_data_dict[L]
        
        # Get the epsilon_log estimates
        epsilon_log_K0_estimate = epsilon_log_K0_estimates[L]
        epsilon_log_K2_estimate = epsilon_log_K2_estimates[L]
        
        # Get mean gamma values from the fixed epsilon_log fits
        gamma_K0_values = [fixed_epsilon_log_K0_fits[L].param[2] for L in L_values]
        mean_gamma_K0 = mean(gamma_K0_values)
        
        gamma_K2_values = [fixed_epsilon_log_K2_fits[L].param[2] for L in L_values]
        mean_gamma_K2 = mean(gamma_K2_values)
        
        # Get fixed epsilon_log fit results
        fixed_fit_K0 = fixed_epsilon_log_K0_fits[L]
        fixed_fit_K2 = fixed_epsilon_log_K2_fits[L]
        
        # Fit K0 data with fixed epsilon_log and gamma
        fixed_epsilon_log_and_gamma_fit_K0 = fit_with_fixed_epsilon_log_and_gamma(
            L, i, fit_K0.param, K0_functional_form, 
            K0_fixed_epsilon_log, K0_fixed_epsilon_log_and_gamma,
            epsilon_log_K0_estimate, mean_gamma_K0, 1, "K=0", addon)
        fixed_epsilon_log_and_gamma_K0_fits[L] = fixed_epsilon_log_and_gamma_fit_K0
        
        # Fit K2 data with fixed epsilon_log and gamma
        fixed_epsilon_log_and_gamma_fit_K2 = fit_with_fixed_epsilon_log_and_gamma(
            L, i, fit_K2.param, K2_functional_form, 
            K2_fixed_epsilon_log, K2_fixed_epsilon_log_and_gamma,
            epsilon_log_K2_estimate, mean_gamma_K2, 4, "K≥2", addon)
        fixed_epsilon_log_and_gamma_K2_fits[L] = fixed_epsilon_log_and_gamma_fit_K2
        
        # Create plots for K0 and K2 fits
        for saddle_type in ["K=0", "K≥2"]
            saddle_index = saddle_type == "K=0" ? 1 : 4
            
            # Get appropriate parameters based on saddle type
            if saddle_type == "K=0"
                fit_params = fit_K0.param
                fixed_fit = fixed_fit_K0
                epsilon_log_estimate = epsilon_log_K0_estimate
                mean_gamma = mean_gamma_K0
                fixed_epsilon_log_and_gamma_fit = fixed_epsilon_log_and_gamma_fit_K0
                fit_func = K0_functional_form
                fixed_epsilon_log_func = K0_fixed_epsilon_log
                fixed_epsilon_log_and_gamma_func = K0_fixed_epsilon_log_and_gamma
            else
                fit_params = fit_K2.param
                fixed_fit = fixed_fit_K2
                epsilon_log_estimate = epsilon_log_K2_estimate
                mean_gamma = mean_gamma_K2
                fixed_epsilon_log_and_gamma_fit = fixed_epsilon_log_and_gamma_fit_K2
                fit_func = K2_functional_form
                fixed_epsilon_log_func = K2_fixed_epsilon_log
                fixed_epsilon_log_and_gamma_func = K2_fixed_epsilon_log_and_gamma
            end
            
            # Create plot
            fixed_epsilon_log_and_gamma_fit_plot = plot(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
                ylabel="$saddle_type Proportions",
                title="L = $L, $saddle_type",
            legend=:topleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05),
            ylims=(0, 1.05)
        )
        
            # Plot the data points
        E0_bin_values, k_saddle_proportions = data_dict[L]
            energy_density_data = -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L)))
            scatter!(fixed_epsilon_log_and_gamma_fit_plot,
                energy_density_data,
                k_saddle_proportions[:, saddle_index],
                label="Data",
                color=:gray,
                marker=:circle,
                markersize=4,
                alpha=0.7
            )
            
            # Plot the original fit
            energy_range = LinRange(-0.7, -0.05, 500)
        plot!(
                fixed_epsilon_log_and_gamma_fit_plot,
                energy_range,
                fit_func(energy_range, fit_params),
                label="Original fit (ϵ_log=$(round(fit_params[2], digits=4)), γ=$(round(fit_params[3], digits=2)))",
            color=:black,
            linewidth=2
        )
        
            # Plot the fixed epsilon_log fit
            plot!(
                fixed_epsilon_log_and_gamma_fit_plot,
                energy_range,
                fixed_epsilon_log_func(energy_range, fixed_fit.param, epsilon_log_estimate),
                label="Fixed ϵ_log fit (ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(fixed_fit.param[2], digits=2)))",
                color=:blue,
                linewidth=2,
                linestyle=:dash
            )

            # Plot the fixed epsilon_log and gamma fit
            plot!(
                fixed_epsilon_log_and_gamma_fit_plot,
                energy_range,
                fixed_epsilon_log_and_gamma_func(energy_range, fixed_epsilon_log_and_gamma_fit.param, epsilon_log_estimate, mean_gamma),
                label="Fixed ϵ_log and γ fit (ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(mean_gamma, digits=2)))",
                color=:red,
                linewidth=2,
                linestyle=:dash
            )
            
            # Add vertical lines for the two epsilon_log values
            vline!(fixed_epsilon_log_and_gamma_fit_plot, [fit_params[2]], linecolor=:black, linestyle=:dash, linewidth=1, 
                   label="Original ϵ_log")
            vline!(fixed_epsilon_log_and_gamma_fit_plot, [epsilon_log_estimate], linecolor=:blue, linestyle=:dash, linewidth=1, 
                   label="Fixed ϵ_log")
            
            # Filter valid data for RMSE calculation
            valid_rows = .!isnan.(k_saddle_proportions[:, saddle_index]) .& 
                        .!isinf.(k_saddle_proportions[:, saddle_index]) .& 
                        .!isnan.(E0_bin_values) .& 
                        .!isinf.(E0_bin_values)
            E0_bin_values_filtered = E0_bin_values[valid_rows]
            k_saddle_proportions_filtered = k_saddle_proportions[valid_rows, saddle_index]
            energy_density = -1.0 .+ (E0_bin_values_filtered./-solved_configuration_energy(RubiksCube(L)))
            
            # Calculate fit quality metrics
            original_rmse = sqrt(sum((fit_func(energy_density, fit_params) .- k_saddle_proportions_filtered).^2) / 
                                length(k_saddle_proportions_filtered))
            fixed_rmse = sqrt(sum((fixed_epsilon_log_func(energy_density, fixed_fit.param, epsilon_log_estimate) .- 
                                    k_saddle_proportions_filtered).^2) / length(k_saddle_proportions_filtered))
            fixed_epsilon_log_and_gamma_rmse = sqrt(sum((fixed_epsilon_log_and_gamma_func(energy_density, fixed_epsilon_log_and_gamma_fit.param, 
                                                        epsilon_log_estimate, mean_gamma) .- k_saddle_proportions_filtered).^2) / 
                                                    length(k_saddle_proportions_filtered))

            # Print fit parameters and quality metrics
            println("L=$L: Original $saddle_type fit: α=$(round(fit_params[1], digits=2)), ϵ_log=$(round(fit_params[2], digits=4)), γ=$(round(fit_params[3], digits=2))")
            println("L=$L: Fixed $saddle_type ϵ_log fit: α=$(round(fixed_fit.param[1], digits=2)), ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(fixed_fit.param[2], digits=2))")
            println("L=$L: Fixed $saddle_type ϵ_log and γ fit: α=$(round(fixed_epsilon_log_and_gamma_fit.param[1], digits=2)), ϵ_log=$(round(epsilon_log_estimate, digits=4)), γ=$(round(mean_gamma, digits=2))")
            println("L=$L: $saddle_type Original RMSE: $(round(original_rmse, digits=5))")
            println("L=$L: $saddle_type Fixed ϵ_log RMSE: $(round(fixed_rmse, digits=5))")
            println("L=$L: $saddle_type Fixed ϵ_log and γ RMSE: $(round(fixed_epsilon_log_and_gamma_rmse, digits=5))")
            
            # Save and display the plot
            saddle_prefix = saddle_type == "K=0" ? "K0" : "K2"
            savefig(fixed_epsilon_log_and_gamma_fit_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(saddle_prefix)_fixed_epsilon_log_and_gamma_fit$(addon).pdf")
            savefig(fixed_epsilon_log_and_gamma_fit_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/L$(L)_$(saddle_prefix)_fixed_epsilon_log_and_gamma_fit$(addon).png")
            display(fixed_epsilon_log_and_gamma_fit_plot)
        end
    end
    
    # Create a 2x2 comparison plot for all parameters with original, fixed epsilon_log, and fixed epsilon_log and gamma fits
    full_params_comparison_plot = plot(
        layout=(2,2),
        size=(900, 800),
        margin=5mm,
    )
    
    # Plot alpha vs L for K0
    plot!(
        full_params_comparison_plot[1],
        [L for L in L_values],
        [fit_data_dict[L][1].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha \textrm{ for } K=0",
        title="Alpha Parameter (K=0)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )
    
    plot!(
        full_params_comparison_plot[1],
        [L for L in L_values],
        [fixed_epsilon_log_K0_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )
    
    plot!(
        full_params_comparison_plot[1],
        [L for L in L_values],
        [fixed_epsilon_log_and_gamma_K0_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log and γ fit",
        marker=:square,
        color=:red
    )
    
    # Plot alpha vs L for K2
    plot!(
        full_params_comparison_plot[2],
        [L for L in L_values],
        [fit_data_dict[L][2].param[1] for L in L_values],
        xlabel=L"L",
        ylabel=L"\alpha \textrm{ for } K \geq 2",
        title="Alpha Parameter (K≥2)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:bottomright
    )

    plot!(
        full_params_comparison_plot[2],
        [L for L in L_values],
        [fixed_epsilon_log_K2_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    plot!(
        full_params_comparison_plot[2],
        [L for L in L_values],
        [fixed_epsilon_log_and_gamma_K2_fits[L].param[1] for L in L_values],
        label="Fixed ϵ_log and γ fit",
        marker=:square,
        color=:red
    )
    
    # Plot gamma vs 1/L for K0
    plot!(
        full_params_comparison_plot[3],
        [1/L for L in L_values],
        [fit_data_dict[L][1].param[3] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\gamma \textrm{ for } K=0",
        title="Gamma Parameter (K=0)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:topright
    )

    plot!(
        full_params_comparison_plot[3],
        [1/L for L in L_values],
        [fixed_epsilon_log_K0_fits[L].param[2] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    mean_gamma_K0 = mean([fixed_epsilon_log_K0_fits[L].param[2] for L in L_values])
    plot!(
        full_params_comparison_plot[3],
        [1/L for L in L_values],
        fill(mean_gamma_K0, length(L_values)),
        label="Fixed ϵ_log and γ fit",
        marker=:square,
        color=:red
    )
    
    # Plot gamma vs 1/L for K2
    plot!(
        full_params_comparison_plot[4],
        [1/L for L in L_values],
        [fit_data_dict[L][2].param[3] for L in L_values],
        xlabel=L"1/L",
        ylabel=L"\gamma \textrm{ for } K \geq 2",
        title="Gamma Parameter (K≥2)",
        label="Original fit",
        marker=:circle,
        color=:black,
        legend=:topright
    )

    plot!(
        full_params_comparison_plot[4],
        [1/L for L in L_values],
        [fixed_epsilon_log_K2_fits[L].param[2] for L in L_values],
        label="Fixed ϵ_log fit",
        marker=:diamond,
        color=:blue
    )

    mean_gamma_K2 = mean([fixed_epsilon_log_K2_fits[L].param[2] for L in L_values])
    plot!(
        full_params_comparison_plot[4],
        [1/L for L in L_values],
        fill(mean_gamma_K2, length(L_values)),
        label="Fixed ϵ_log and γ fit",
        marker=:square,
        color=:red
    )
    
    # Fit lines for alpha parameters with all three types of fits
    L_vals = [L for L in L_values]
    
    # K0 fits
    alpha_K0_orig_vals = [fit_data_dict[L][1].param[1] for L in L_values]
    fit_alpha_K0_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_orig_vals, [0.0, 0.0])
    
    alpha_K0_fixed_vals = [fixed_epsilon_log_K0_fits[L].param[1] for L in L_values]
    fit_alpha_K0_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_fixed_vals, [0.0, 0.0])
    
    alpha_K0_fixed_gamma_vals = [fixed_epsilon_log_and_gamma_K0_fits[L].param[1] for L in L_values]
    fit_alpha_K0_fixed_gamma = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K0_fixed_gamma_vals, [0.0, 0.0])
    
    # K2 fits
    alpha_K2_orig_vals = [fit_data_dict[L][2].param[1] for L in L_values]
    fit_alpha_K2_orig = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_orig_vals, [0.0, 0.0])
    
    alpha_K2_fixed_vals = [fixed_epsilon_log_K2_fits[L].param[1] for L in L_values]
    fit_alpha_K2_fixed = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_fixed_vals, [0.0, 0.0])
    
    alpha_K2_fixed_gamma_vals = [fixed_epsilon_log_and_gamma_K2_fits[L].param[1] for L in L_values]
    fit_alpha_K2_fixed_gamma = curve_fit((x, p) -> p[1] .+ p[2].*x, L_vals, alpha_K2_fixed_gamma_vals, [0.0, 0.0])
    
    # Plot fit lines
    x_fit_L = LinRange(min(L_values...), max(L_values...), 100)
    
    # K0 alpha fit lines
    plot!(
        full_params_comparison_plot[1],
        x_fit_L,
        fit_alpha_K0_orig.param[1] .+ fit_alpha_K0_orig.param[2].*x_fit_L,
        label="Original: $(round(fit_alpha_K0_orig.param[1], digits=2)) + $(round(fit_alpha_K0_orig.param[2], digits=2))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        full_params_comparison_plot[1],
        x_fit_L,
        fit_alpha_K0_fixed.param[1] .+ fit_alpha_K0_fixed.param[2].*x_fit_L,
        label="Fixed ϵ_log: $(round(fit_alpha_K0_fixed.param[1], digits=2)) + $(round(fit_alpha_K0_fixed.param[2], digits=2))·L",
        color=:blue,
        linestyle=:dash
    )
    
    plot!(
        full_params_comparison_plot[1],
        x_fit_L,
        fit_alpha_K0_fixed_gamma.param[1] .+ fit_alpha_K0_fixed_gamma.param[2].*x_fit_L,
        label="Fixed ϵ_log & γ: $(round(fit_alpha_K0_fixed_gamma.param[1], digits=2)) + $(round(fit_alpha_K0_fixed_gamma.param[2], digits=2))·L",
        color=:red,
        linestyle=:dash
    )
    
    # K2 alpha fit lines
    plot!(
        full_params_comparison_plot[2],
        x_fit_L,
        fit_alpha_K2_orig.param[1] .+ fit_alpha_K2_orig.param[2].*x_fit_L,
        label="Original: $(round(fit_alpha_K2_orig.param[1], digits=2)) + $(round(fit_alpha_K2_orig.param[2], digits=2))·L",
        color=:black,
        linestyle=:dash
    )

    plot!(
        full_params_comparison_plot[2],
        x_fit_L,
        fit_alpha_K2_fixed.param[1] .+ fit_alpha_K2_fixed.param[2].*x_fit_L,
        label="Fixed ϵ_log: $(round(fit_alpha_K2_fixed.param[1], digits=2)) + $(round(fit_alpha_K2_fixed.param[2], digits=2))·L",
        color=:blue,
        linestyle=:dash
    )

    plot!(
        full_params_comparison_plot[2],
        x_fit_L,
        fit_alpha_K2_fixed_gamma.param[1] .+ fit_alpha_K2_fixed_gamma.param[2].*x_fit_L,
        label="Fixed ϵ_log & γ: $(round(fit_alpha_K2_fixed_gamma.param[1], digits=2)) + $(round(fit_alpha_K2_fixed_gamma.param[2], digits=2))·L",
        color=:red,
        linestyle=:dash
    )

    # Print fit results
    println("K0 alpha fits:")
    println("  Original: $(round(fit_alpha_K0_orig.param[1], digits=4)) + $(round(fit_alpha_K0_orig.param[2], digits=4))·L")
    println("  Fixed ϵ_log: $(round(fit_alpha_K0_fixed.param[1], digits=4)) + $(round(fit_alpha_K0_fixed.param[2], digits=4))·L")
    println("  Fixed ϵ_log & γ: $(round(fit_alpha_K0_fixed_gamma.param[1], digits=4)) + $(round(fit_alpha_K0_fixed_gamma.param[2], digits=4))·L")
    
    println("K2 alpha fits:")
    println("  Original: $(round(fit_alpha_K2_orig.param[1], digits=4)) + $(round(fit_alpha_K2_orig.param[2], digits=4))·L")
    println("  Fixed ϵ_log: $(round(fit_alpha_K2_fixed.param[1], digits=4)) + $(round(fit_alpha_K2_fixed.param[2], digits=4))·L")
    println("  Fixed ϵ_log & γ: $(round(fit_alpha_K2_fixed_gamma.param[1], digits=4)) + $(round(fit_alpha_K2_fixed_gamma.param[2], digits=4))·L")
    
    # Save and display the comparison plot
    savefig(full_params_comparison_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_all_fixed_vs_original$(addon).pdf")
    savefig(full_params_comparison_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/parameters_all_fixed_vs_original$(addon).png")
    display(full_params_comparison_plot)
    
    
    ##### COMBINED OVERLAY PLOTS FOR K=0 AND K>=2 FITS #####
    
    # Create combined plot with K=0 and K>=2 fixed epsilon_log fits
    combined_fixed_epsilon_log_plot = plot(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Saddle Index Proportions",
        legend=:bottomleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05),
        ylims=(0, 1.05)
    )
    
    # Create combined plot with K=0 and K>=2 fixed epsilon_log and gamma fits
    combined_fixed_epsilon_log_and_gamma_plot = plot(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Saddle Index Proportions",
        legend=:bottomleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05),
        ylims=(0, 1.05)
    )
    
    energy_range = LinRange(-0.7, -0.05, 500)
    
    for (i, L) in enumerate(L_values)
        # Get the data for this L value
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Get energy density data
        energy_density_data = -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L)))
        
        # Get the fixed parameters
        epsilon_log_K0_estimate = epsilon_log_K0_estimates[L]
        epsilon_log_K2_estimate = epsilon_log_K2_estimates[L]
        
        mean_gamma_K0 = mean([fixed_epsilon_log_K0_fits[L].param[2] for L in L_values])
        mean_gamma_K2 = mean([fixed_epsilon_log_K2_fits[L].param[2] for L in L_values])
        
        # Plot K=0 data points
        scatter!(combined_fixed_epsilon_log_plot,
            energy_density_data,
            k_saddle_proportions[:, 1],
            label=i==1 ? "K=0" : "",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot,
            energy_density_data,
            k_saddle_proportions[:, 1],
            label=i==1 ? "K=0" : "",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
        
        # Plot K≥2 data points
        scatter!(combined_fixed_epsilon_log_plot,
            energy_density_data,
            k_saddle_proportions[:, 4],
            label=i==1 ? "K≥2" : "",
            color=colors[i],
            marker=:diamond,
            markersize=4
        )
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot,
            energy_density_data,
            k_saddle_proportions[:, 4],
            label=i==1 ? "K≥2" : "",
            color=colors[i],
            marker=:diamond,
            markersize=4
        )
        
        # Plot K=0 fixed epsilon_log fit
        plot!(
            combined_fixed_epsilon_log_plot,
            energy_range,
            K0_fixed_epsilon_log(energy_range, fixed_epsilon_log_K0_fits[L].param, epsilon_log_K0_estimate),
            label="",
            color=colors[i],
            linewidth=2,
            linestyle=:dash
        )
        
        # Plot K≥2 fixed epsilon_log fit
        plot!(
            combined_fixed_epsilon_log_plot,
            energy_range,
            K2_fixed_epsilon_log(energy_range, fixed_epsilon_log_K2_fits[L].param, epsilon_log_K2_estimate),
            label="",
            color=colors[i],
            linewidth=2,
            linestyle=:dot
        )
        
        # Plot K=0, fixed epsilon_log and gamma fit
        plot!(
            combined_fixed_epsilon_log_and_gamma_plot,
            energy_range,
            K0_fixed_epsilon_log_and_gamma(energy_range, fixed_epsilon_log_and_gamma_K0_fits[L].param, epsilon_log_K0_estimate, mean_gamma_K0),
            label="",
            color=colors[i],
            linewidth=2,
            linestyle=:dash
        )
        
        # Plot K≥2, fixed epsilon_log and gamma fit
        plot!(
            combined_fixed_epsilon_log_and_gamma_plot,
            energy_range,
            K2_fixed_epsilon_log_and_gamma(energy_range, fixed_epsilon_log_and_gamma_K2_fits[L].param, epsilon_log_K2_estimate, mean_gamma_K2),
            label="",
            color=colors[i],
            linewidth=2,
            linestyle=:dot
        )
    end
    
    # Add L value entries to legend for both plots
    for (i, L) in enumerate(L_values)
        scatter!(combined_fixed_epsilon_log_plot, [], [], 
                label="L = $L", color=colors[i], marker=:square, markersize=4)
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot, [], [], 
                label="L = $L", color=colors[i], marker=:square, markersize=4)
    end
    
    # Add line style indicators to legend
    plot!(combined_fixed_epsilon_log_plot, [], [], 
          color=:black, linestyle=:dash, label="K=0 fit", linewidth=2)
    plot!(combined_fixed_epsilon_log_plot, [], [], 
          color=:black, linestyle=:dot, label="K≥2 fit", linewidth=2)
    
    plot!(combined_fixed_epsilon_log_and_gamma_plot, [], [], 
          color=:black, linestyle=:dash, label="K=0 fit", linewidth=2)
    plot!(combined_fixed_epsilon_log_and_gamma_plot, [], [], 
          color=:black, linestyle=:dot, label="K≥2 fit", linewidth=2)
    
    # Save and display the combined plots
    savefig(combined_fixed_epsilon_log_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K0_K2_fixed_epsilon_log_fits$(addon).pdf")
    savefig(combined_fixed_epsilon_log_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K0_K2_fixed_epsilon_log_fits$(addon).png")
    display(combined_fixed_epsilon_log_plot)
    
    savefig(combined_fixed_epsilon_log_and_gamma_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K0_K2_fixed_epsilon_log_and_gamma_fits$(addon).pdf")
    savefig(combined_fixed_epsilon_log_and_gamma_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_K0_K2_fixed_epsilon_log_and_gamma_fits$(addon).png")
    display(combined_fixed_epsilon_log_and_gamma_plot)
    

    ##### #####













    ##### FIXED EPSILON LOG AND GAMMA COLLAPSED PLOT #####
    # Make a combined collapsed plot of the K>=2 saddle index proportions but plotted against (a_0 + a_1L)(epsilon - epsilon_log)
    # Where a_0 and a_1 are the parameters of the fixed epsilon_log and gamma fit of alpha

    # Get the parameters of the fixed epsilon_log and gamma fit of alpha for K2
    a_0 = fit_alpha_K2_fixed_gamma.param[1]
    a_1 = fit_alpha_K2_fixed_gamma.param[2]




    collapsed_plot = plot(
        xlabel="($(round(a_0, digits=2)) + $(round(a_1, digits=2))L)"*"(ϵ - ϵ_log(L))",
        ylabel="K≥2 Proportions",
        title="Collapsed Plot",
        legend=:topleft,
    )


    for (i, L) in enumerate(L_values)

        E0_bin_values, k_saddle_proportions = data_dict[L]

        energy_density_data = -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L)))

        # Filter energy densities that are within -1 and 0
        valid_rows = energy_density_data .>= -1.0 .&& energy_density_data .<= 0.0
        energy_density_data_filtered = energy_density_data[valid_rows]
        k_saddle_proportions_filtered = k_saddle_proportions[valid_rows, 4]

        # Get the epsilon_log values for K2
        epsilon_log_value = epsilon_log_K2_estimates[L]

        # Get the collapsed values
        collapsed_values = (a_0 .+ a_1.*L).*(energy_density_data_filtered .- epsilon_log_value)

        scatter!(
            collapsed_plot,
            collapsed_values,
            k_saddle_proportions_filtered,
            label="L = $L",
            color=colors[i],
        )
    end
    
    
    savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_fixed_epsilon_log_and_gamma$(addon).pdf")
    savefig(collapsed_plot, "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/collapsed_fixed_epsilon_log_and_gamma$(addon).png")
    display(collapsed_plot)
    

    ##### #####









    ##### COMBINED PLOTS FOR FIXED FITS #####
    
    # Create a combined plot for all epsilon_log fixed fits
    combined_fixed_epsilon_log_plot = plot(
            xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K = 0"*" and "*L"K \geq 2"*" Proportions",
        legend=:bottomleft,
            xguidefontsize=12,
            yguidefontsize=12,
            margin=5mm,
            xlims=(-0.7, -0.05),
            ylims=(0, 1.05)
        )
        
    
    # Create a combined plot for all epsilon_log and gamma fixed fits
    combined_fixed_epsilon_log_and_gamma_plot = plot(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K = 0"*" and "*L"K \geq 2"*" Proportions",
        legend=:bottomleft,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05),
        ylims=(0, 1.05)
    )
    
    # Get mean gamma values from all fixed epsilon_log fits for K0 and K2
    gamma_K0_values = [fixed_epsilon_log_K0_fits[L].param[2] for L in L_values]
    mean_gamma_K0 = mean(gamma_K0_values)
    
    gamma_K2_values = [fixed_epsilon_log_K2_fits[L].param[2] for L in L_values]
    mean_gamma_K2 = mean(gamma_K2_values)
    
        energy_range = LinRange(-0.7, -0.05, 500)
    
    for (i, L) in enumerate(L_values)
        # Get the data for this L value
        E0_bin_values, k_saddle_proportions = data_dict[L]
        
        # Transform to energy density
        energy_density_data = -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L)))
        
        # Get the fixed parameters for epsilon_log (using K2 estimates for consistency)
        epsilon_log_estimate = epsilon_log_K2_estimates[L]
        
        # Plot K=0 data points
        scatter!(combined_fixed_epsilon_log_plot,
            energy_density_data,
            k_saddle_proportions[:, 1],
            label=i==1 ? "K=0" : "",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot,
            energy_density_data,
            k_saddle_proportions[:, 1],
            label=i==1 ? "K=0" : "",
            color=colors[i],
            marker=:circle,
            markersize=4
        )
        
        # Plot K≥2 data points
        scatter!(combined_fixed_epsilon_log_plot,
            energy_density_data,
            k_saddle_proportions[:, 4],
            label=i==1 ? "K≥2" : "",
            color=colors[i],
            marker=:diamond,
            markersize=4
        )
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot,
            energy_density_data,
            k_saddle_proportions[:, 4],
            label=i==1 ? "K≥2" : "",
            color=colors[i],
            marker=:diamond,
            markersize=4
        )
        
        # Plot K≥2 fixed epsilon_log fit only
        plot!(
            combined_fixed_epsilon_log_plot,
            energy_range,
            K2_fixed_epsilon_log(energy_range, fixed_epsilon_log_K2_fits[L].param, epsilon_log_estimate),
            label=i==1 ? "K≥2, fit" : "",
            color=colors[i],
            linewidth=2,
            linestyle=:dash
        )

        # Plot K≥2 fixed epsilon_log and gamma fit only
        plot!(
            combined_fixed_epsilon_log_and_gamma_plot,
            energy_range,
            K2_fixed_epsilon_log_and_gamma(energy_range, fixed_epsilon_log_and_gamma_K2_fits[L].param, epsilon_log_estimate, mean_gamma_K2),
            label=i==1 ? "K≥2, fit" : "",
            color=colors[i],
            linewidth=2,
            linestyle=:dash
        )
    end
    
    
    # Add L value entries to legend
    for (i, L) in enumerate(L_values)
        scatter!(combined_fixed_epsilon_log_plot, [], [], 
                label="L = $L", color=colors[i], marker=:square, markersize=4)
        
        scatter!(combined_fixed_epsilon_log_and_gamma_plot, [], [], 
                label="L = $L", color=colors[i], marker=:square, markersize=4)
    end
    
    # Save and display the combined plots
    savefig(combined_fixed_epsilon_log_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_fixed_epsilon_log_fits$(addon).pdf")
    savefig(combined_fixed_epsilon_log_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_fixed_epsilon_log_fits$(addon).png")
    display(combined_fixed_epsilon_log_plot)
    
    savefig(combined_fixed_epsilon_log_and_gamma_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_fixed_epsilon_log_and_gamma_fits$(addon).pdf")
    savefig(combined_fixed_epsilon_log_and_gamma_plot, 
            "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_fixed_epsilon_log_and_gamma_fits$(addon).png")
    display(combined_fixed_epsilon_log_and_gamma_plot)
end
































# all_saddle_index_proportions_figure_by_E1_E0_L_analysis([13,11,9,7,5,3])


# combined_saddle_index_proportions_figure_by_L([13,11,9,7,5,3]; include_K1_in_K0=true)
combined_saddle_index_proportions_figure_by_L([15,13,11,9,7,5,3])
































