import Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

using LaTeXStrings
using DelimitedFiles
using Plots
using Plots.PlotMeasures
using Colors
using Statistics
using LsqFit
using Distributions

include("../core/rubiks_cube.jl")

linear_model(x, p) = p[1] .* x .+ p[2]

function relaxed_anneal_figure(extraction::Bool=false, annotations=false)

    output_name = "L=11_relaxed_anneal_figure_raw"

    # --- L=11 Figure ---
    models = ["clean", "inherent_disorder"]
    Ls = [11]
    swap_move_probabilities = [0.0, 1.0]
    trials = 50



    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)


    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                if L==11 && swap_move_probability == 1.0 && model == "inherent_disorder"
                    N_T = 100
                else
                    N_T = 200
                end


                temperatures = zeros(N_T)
                running_total_energy_densities_by_temperature = zeros(N_T)
                running_total_squared_energies_by_temperature = zeros(N_T)


                actual_number_of_trials=0
                for trial in 1:trials
                    filename = "results/relaxed_anneal_results/data/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_energy_densities_by_temperature .+= data_matrix[:,3]
                        running_total_squared_energies_by_temperature .+= data_matrix[:,4]
                        
                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_energy_densities_by_temperature/actual_number_of_trials
                
                average_squared_energies_by_temperature = running_total_squared_energies_by_temperature / actual_number_of_trials
                normalization_factor = 12 * L * (L - 1)
                average_squared_energy_densities_by_temperature = average_squared_energies_by_temperature ./ (normalization_factor^2)

                # Calculate standard errors using variance: SE = sqrt(Var(E)/(N_trials * N_average))
                # Where Var(E) = <E²> - <E>²
                N_average_per_temperature = 100
                variances = average_squared_energy_densities_by_temperature .- results_dictionary[(model, L, swap_move_probability, "average_energy_densities")].^2
                standard_errors = sqrt.(variances ./ (actual_number_of_trials * N_average_per_temperature))
                results_dictionary[(model, L, swap_move_probability, "standard_errors")] = standard_errors            
            end
        end
    end


    ### --- PLOT DATA ---
    graph = plot(title="", xlabel="Temperature, "*L"T", legend=:bottomright, yaxis="Average Energy Density, "*L"\langle\! \epsilon \rangle = \langle\! E/|\!\!E_s|\!\rangle", ylims=(-1.0,-0.1))


    color_index = 1
    for model in models
        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_energy_densities = results_dictionary[(model, L, swap_move_probability, "average_energy_densities")]
                standard_errors = results_dictionary[(model, L, swap_move_probability, "standard_errors")]

                
                marker_style = swap_move_probability == 0.0 ? :diamond : :circle
                marker_size = 3

                label = ""
                if swap_move_probability == 1.0
                    label *= "L = $L"
                    if model == "inherent_disorder" && swap_move_probability == 1.0
                        label *= ", Randomised"
                    end
                    if model == "custom" && swap_move_probability == 1.0
                        label *= ", Randomisation Instance"
                    end
                end

                # Plot with data points and error bars 
                scatter!(graph, temperatures, average_energy_densities, 
                       yerror=standard_errors,
                       label=label, color=colors[mod1(color_index,length(colors))], 
                       markersize=marker_size, markershape=marker_style, markerstrokewidth=0.05)
            end
            color_index += 1
        end
    end

    # ### --- ADD ANNOTATIONS TO GRAPH ---
    # if annotations
    #     if "inherent_disorder" in models
    #         # annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"\bar{\epsilon}^*", 12, alex_red, ))])

    #         T_star = 0.8
    #         # T_on = 2.25

    #         println("T* Inherent Disorder = $T_star")
    #         # println("T_on Inherent Disorder = $T_on")
    #         println("epsilon at T* Inherent Disorder = $(results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_star))])")
    #         # annotate!(graph, [(0.6, ylims(graph)[2]-0.18, Plots.text(L"(\bar{T}^{*}\!\!\!\!,\bar{\epsilon}^{\!\!\!*}\!\!)", 12, alex_red))])
    #         # annotate!(graph, [(T_on + 0.2, ylims(graph)[2]-0.06, Plots.text(L"(\bar{T}^{on}\!\!\!\!,\bar{\epsilon}^{on})", 12, alex_red))])


    #         # Add dot at (T_star, inherent_disorder average_energy_densities at T nearest to T_star)
    #         inherent_disorder_average_energy_densities = results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")]
    #         scatter!([T_star], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_star))]], color=alex_red, markerstrokecolor=alex_red, label="", markersize=2.5)
                    

    #         # Get \epsilon_on as the average energy density  at T_on
    #         # inherent_disorder_epsilon_on = inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_on))]
    #         # println("Inherent Disorder epsilon_on = $inherent_disorder_epsilon_on")        
    #         # Add dot at (T_on, inherent_disorder average_energy_densities at T nearest to T_on)
    #         # scatter!([T_on], [inherent_disorder_epsilon_on], color=alex_red, label="", markersize=1.5)

    #         if extraction

    #             # Add dot at (lowest T value with average_energy_densities values, inherent_disorder average_energy_densities at T nearest to 0)
    #             # scatter!([minimum(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")])], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")]))]], color=alex_red, label="", markersize=1.5)

    #             inherent_disorder_epsilon_star_avg = minimum(results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")])
    #             inherent_disorder_epsilon_0_avg = minimum(results_dictionary[("inherent_disorder", 11, 1.0, "average_energy_densities")])

    #             println("Inherent Disorder epsilon^* = $inherent_disorder_epsilon_star_avg")
    #             println("Inherent Disorder epsilon_0 = $inherent_disorder_epsilon_0_avg")



    #         end
    #     end

    #     if "clean" in models
    #         if extraction
    #             annotate!(graph, [(0.3, ylims(graph)[2]-0.35, Plots.text(L"{\epsilon}^*", 12, alex_alt_blue, ))])
    #         end

    #         T_c = 0.9
    #         println("T_c = $T_c")
    #         # annotate!(graph, [(T_c+0.25, ylims(graph)[1]+0.05, Plots.text(L"T_c", 12, alex_alt_blue))])

    #         T_star = 0.92
    #         println("T* Clean = $T_star")
    #         println("epsilon at T* Clean = $(results_dictionary[("clean", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")].-T_star))])")
    #         if extraction
    #             annotate!(graph, [(T_star+0.6, ylims(graph)[2]-0.3, Plots.text(L"T^*", 12, alex_alt_blue))])
    #         end

    #         if extraction
    #             # Add dot at (T_star, clean_average_energy_densities at T nearest to T_star)
    #             clean_average_energy_densities = results_dictionary[("clean", 11, 0.0, "average_energy_densities")]
    #             scatter!([T_star], [clean_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")].-T_star))]], color=alex_alt_blue, label="", markersize=1.5)
            
    #             # Add dot at (lowest T value with average_energy_densities values, clean_average_energy_densities at T nearest to 0)
    #             scatter!([minimum(results_dictionary[("clean", 11, 0.0, "temperatures")])], [clean_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")]))]], color=alex_alt_blue, label="", markersize=1.5)

    #             # Add dot at (T_c, clean-average_energy_densities at T nearest to T_c)
    #             clean_swap_average_energy_densities = results_dictionary[("clean", 11, 1.0, "average_energy_densities")]
    #             scatter!([T_c], [clean_swap_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 1.0, "temperatures")].-T_c))]], color=alex_alt_blue, label="", markersize=1.5)

    #             clean_epsilon_star_avg = minimum(results_dictionary[("clean", 11, 0.0, "average_energy_densities")])
    #             clean_epsilon_0_avg = minimum(results_dictionary[("clean", 11, 1.0, "average_energy_densities")])

    #             println("Clean epsilon^* = $clean_epsilon_star_avg")
    #             println("Clean epsilon_0 = $clean_epsilon_0_avg")
    #         end
    #     end

    #     if "custom" in models
    #         annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"{\epsilon}^*", 12, alex_blue, ))])

    #         T_star = 0.86
    #         println("T* Custom = $T_star")
    #         println("epsilon at T* Custom = $(results_dictionary[("custom", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))])")
    #         annotate!(graph, [(T_star+0.4, ylims(graph)[2]-0.28, Plots.text(L"T^*", 12, alex_blue))])
    #         if extraction
    #             # Add dot at (T_star, custom average_energy_densities at T nearest to T_star)
    #             custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
    #             scatter!([T_star], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))]], color=alex_blue, label="", markersize=1.5)

    #             # Add dot at (lowest T value with average_energy_densities values, custom average_energy_densities at T nearest to 0)
    #             custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
    #             scatter!([minimum(results_dictionary[("custom", 11, 0.0, "temperatures")])], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")]))]], color=alex_blue, label="", markersize=1.5)
            
                        
    #             custom_epsilon_star_avg = minimum(results_dictionary[("custom", 11, 0.0, "average_energy_densities")])
    #             custom_epsilon_0_avg = minimum(results_dictionary[("custom", 11, 1.0, "average_energy_densities")])

    #             println("Custom epsilon^* = $custom_epsilon_star_avg")
    #             println("Custom epsilon_0 = $custom_epsilon_0_avg")
    #         end
    #     end
    # end



 

    ### --- SAVE GRAPH ---
    # savefig(graph, "results/final_paper_results/$(output_name).png")
    # savefig(graph, "results/final_paper_results/$(output_name).pdf")
    savefig(graph, "results/relaxed_anneal_results/$(output_name).png")
    savefig(graph, "results/relaxed_anneal_results/$(output_name).pdf")
    display(graph)

    ### --- PRINT ERRORS ---
    println("The following files do not exist:")
    for filename in filenames_that_do_not_exist
        println(filename)
    end

end

























function all_L_relaxed_anneal_figure(model::String="inherent_disorder")

    output_name = "all_L_relaxed_anneal_figure_$(model)"

    ## --- All L Figure ---
    models = [model]
    Ls = [5, 7, 9, 11]
    swap_move_probabilities = [0.0, 1.0]
    trials = 50
    N_T = 200




    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)


    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                if L==11 && swap_move_probability == 1.0 && model == "inherent_disorder"
                    N_T = 100
                else
                    N_T = 200
                end


                temperatures = zeros(N_T)
                running_total_energy_densities_by_temperature = zeros(N_T)
                running_total_squared_energies_by_temperature = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    filename = "results/relaxed_anneal_results/data/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_energy_densities_by_temperature .+= data_matrix[:,3]
                        running_total_squared_energies_by_temperature .+= data_matrix[:,4]
                        
                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_energy_densities_by_temperature/actual_number_of_trials
                
                average_squared_energies_by_temperature = running_total_squared_energies_by_temperature / actual_number_of_trials
                normalization_factor = 12 * L * (L - 1)
                average_squared_energy_densities_by_temperature = average_squared_energies_by_temperature ./ (normalization_factor^2)
                
                # Calculate standard errors using variance: SE = sqrt(Var(E)/(N_trials * N_average))
                # Where Var(E) = <E²> - <E>²
                N_average_per_temperature = 100
                variances = average_squared_energy_densities_by_temperature .- results_dictionary[(model, L, swap_move_probability, "average_energy_densities")].^2
                standard_errors = sqrt.(variances ./ (actual_number_of_trials * N_average_per_temperature))
                results_dictionary[(model, L, swap_move_probability, "standard_errors")] = standard_errors
            end
        end
    end






    ### --- PLOT DATA ---
    graph = plot(title="", xlabel="Temperature, "*L"T", legend=:bottomright, yaxis="Average Energy Density, "*L"\langle\! \epsilon \rangle = \langle\! E/|\!\!E_s|\!\rangle", ylims=(-1.0,-0.1))

    epsilon_parameters_graph = plot(title="", xlabel=L"1/L", legend=:bottomleft, yaxis="Energy Density, "*L"\epsilon", ylims=(-0.5,-0.19))

    color_index = 1
    for model in models
        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_energy_densities = results_dictionary[(model, L, swap_move_probability, "average_energy_densities")]
                standard_errors = results_dictionary[(model, L, swap_move_probability, "standard_errors")]

                linestyle = swap_move_probability == 0.0 ? :dash : :solid
                marker_style = swap_move_probability == 0.0 ? :diamond : :circle
                marker_size = 3

                label = ""
                if swap_move_probability == 1.0
                    label *= "L = $L"
                    if model == "inherent_disorder" && swap_move_probability == 1.0
                        label *= ", Randomised"
                    end
                    if model == "custom" && swap_move_probability == 1.0
                        label *= ", Randomisation Instance"
                    end
                end

                # TODO possibly redo this section
                # Plot line for the average values
                plot!(graph, temperatures, average_energy_densities, 
                     label=label, color=colors[mod1(color_index,length(colors))], 
                     linestyle=linestyle)
                
                # TODO possibly remove
                # # Add scatter points with error bars at a sparser interval to avoid overcrowding
                # interval = 1 # Only show error bars every 10 points
                # scatter_indices = 1:interval:length(temperatures)
                # scatter!(graph, temperatures[scatter_indices], average_energy_densities[scatter_indices],
                #        yerror=standard_errors[scatter_indices],
                #        label="", color=colors[mod1(color_index,length(colors))], 
                #        markersize=marker_size, markershape=marker_style)

                # Also print ground state energy from minimum energy from ps=1.0
                if swap_move_probability == 1.0
                    println("L = $L, model = $model, ground state energy = $(minimum(average_energy_densities))\n")
                end
            end
            color_index += 1
        end
    end





    ### --- ADD INSET WITH \epsilon^* for all L ---
    # Epsilon star is energy density at T=0

    # epsilon_star_values_clean = Float64[]
    # epsilon_star_values_inherent_disorder = Float64[]
    epsilon_star_values = Float64[]

    for L in Ls

        # average_energy_densities_clean = results_dictionary[("clean", L, 0.0, "average_energy_densities")]
        # average_energy_densities_inherent_disorder = results_dictionary[("inherent_disorder", L, 0.0, "average_energy_densities")]
        average_energy_densities = results_dictionary[(model, L, 0.0, "average_energy_densities")]

        # epsilon_star_clean = minimum(average_energy_densities_clean)
        # epsilon_star_inherent_disorder = minimum(average_energy_densities_inherent_disorder)
        epsilon_star = minimum(average_energy_densities)
        
        # push!(epsilon_star_values_clean, epsilon_star_clean)
        # push!(epsilon_star_values_inherent_disorder, epsilon_star_inherent_disorder)
        push!(epsilon_star_values, epsilon_star)

    end


    ### --- ADD INSET WITH epsilon^* for all L ---

    ylims = model == "inherent_disorder" ? (-0.5,-0.19) : (-0.65,-0.2)
    scatter!(graph, [1/L for L in Ls], epsilon_star_values, inset=bbox(0.3,0.43,0.35,0.35), subplot=2, xlabel=L"1/L", ylabel=L"\epsilon^*", yguidefontsize=12,xguidefontsize=12, color=:green, label=L"\epsilon^*", legend=:bottomleft, xticks=[0.0, 0.05, 0.1,0.15,0.2], ylims=ylims) # Clean (-0.62,-0.4) 
    # Add linear fit to inset for epsilon_star
    fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_star_values, [1.0, 0.0])
    m, c = fit.param
    confidence_intervals = confidence_interval(fit, 0.05)     # Calculate the 95% confidence intervals for the fit parameters
    c_ci = confidence_intervals[2]
    c_error = c - c_ci[1]
    fit_1_L_values = LinRange(0.0, 0.2, 100)
    plot!(graph, fit_1_L_values, m*fit_1_L_values .+ c, subplot=2, color=:green, label="", linestyle=:dash, legend_font=font(6))
    scatter!(graph, [0.0], [c], yerror=(c_error), subplot=2, color=:green, markerstrokecolor=:green, markersize=1, label="", errorbar_color=:green)
    println("m = $m, c = $c for epsilon_star")



    ### --- ADD INSET WITH T^* for all L ---

    # Taken by eye from start of plateau in graph
    T_star_values_clean_dict = Dict(3 => 0.67, 5 => 0.67, 7 => 0.72, 9 => 0.75, 11 => 0.92)
    T_star_values_inherent_disorder_dict = Dict(3 => 0.3, 5 => 0.57, 7 => 0.65, 9 => 0.68, 11 => 0.8)

    T_star_dict = model == "inherent_disorder" ? T_star_values_inherent_disorder_dict : T_star_values_clean_dict
    T_star_values = [T_star_dict[L] for L in Ls]

    ylims = model == "inherent_disorder" ? (0.55,0.95) : (0.65, 1.05) # (0.55,0.82) : (0.65,0.95)
    scatter!(graph, [1/L for L in Ls], T_star_values, inset=bbox(0.67,0.15,0.28,0.28), subplot=3, xlabel=L"1/L", ylabel=L"T^*", yguidefontsize=12,xguidefontsize=12, color=alex_alt_blue, label="", legend=:bottomright, xticks=[0.0, 0.05, 0.1,0.15,0.2], ylims=ylims)
    
    # Add linear fit to inset for T_star
    fit = curve_fit(linear_model, [1/L for L in Ls], T_star_values, [1.0, 0.0])
    m, c = fit.param
    confidence_intervals = confidence_interval(fit, 0.05)
    c_ci = confidence_intervals[2]
    c_error = c - c_ci[1]
    fit_1_L_values = LinRange(0.0, 0.2, 100)
    plot!(graph, fit_1_L_values, m*fit_1_L_values .+ c, subplot=3, color=alex_alt_blue, label="", linestyle=:dash)
    # scatter!(graph, [0.0], [c], yerror=(c_error), subplot=3, color=alex_green, markerstrokecolor=alex_green, markersize=1, label="", errorbar_color=alex_green)
    println("m = $m, c = $c for T_star")











    # CREATE EPSILON PARAMETERS GRAPH IF MODEL IS INHERENT DISORDER
    if model == "inherent_disorder"


        # epsilon^* values
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], epsilon_star_values, color=:green, label=L"\bar\epsilon^*", legend=:bottomleft, xticks=[0.0, 0.05, 0.1,0.15,0.2])
        fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_star_values, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=:green, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error), color=:green, markerstrokecolor=:green, markersize=1, label="", errorbar_color=:green)







        # epsilon_times values
        epsilon_times_dict = Dict(5 => -0.48163802814018375, 7 => -0.4408591687142211, 11 => -0.38709136614808787, 9 => -0.4109476667396038, 3 => -0.5079313805267485)
        epsilon_times_values = [epsilon_times_dict[L] for L in Ls]
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], epsilon_times_values, color=alex_blue, label=L"\bar\epsilon_{\times}", markersize=4)
        fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_times_values, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=alex_blue, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error), color=alex_blue, markerstrokecolor=alex_blue, markersize=1, label="", errorbar_color=alex_blue)

        
        

        # epsilon_int values
        intersection_points_dict = Dict(5 => (-0.3136219162455117, 0.35311899934539354), 7 => (-0.3092084591537035, 0.345445539802899), 11 => (-0.29372816498023085, 0.3461351961459981), 9 => (-0.30103114454968694, 0.3474276292535396), 3 => (-0.3047545163510452, 0.3029234189450653))
        intersection_points = [intersection_points_dict[L][1] for L in Ls]
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], intersection_points, color=:black, label=L"\bar\epsilon_{\rm int}", markersize=4)
        fit = curve_fit(linear_model, [1/L for L in Ls], intersection_points, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=:black, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error), color=:black, markerstrokecolor=:black, markersize=1, label="", errorbar_color=:black)

        


        # Add epsilon_log values to inset
        epsilon_log_dict = Dict(11 => -0.22642192976237527, 9 => -0.21804406645607335, 7 => -0.21403618914215625, 5 => -0.2041995980294149)
        epsilon_log_values = [epsilon_log_dict[L] for L in Ls]
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], epsilon_log_values, color=alex_green, label=L"\bar\epsilon_{\log}", markersize=4)
        fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_log_values, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=alex_green, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error), color=alex_green, markerstrokecolor=alex_green, markersize=1, label="", errorbar_color=alex_green)





        # Add epsilon_otimes values 
        epsilon_otimes_dict = Dict(5 => -0.2137176988190482, 7 => -0.23134348655947637, 11 => -0.23927293015680562, 9 => -0.23676323523181977, 3 => -0.16330243097243677)
        epsilon_otimes_values = [epsilon_otimes_dict[L] for L in Ls]
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], epsilon_otimes_values, color=alex_orange, label=L"\bar\epsilon_{\otimes}", markersize=4)
        fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_otimes_values, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=alex_orange, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error), color=alex_orange, markerstrokecolor=alex_orange, markersize=1, label="", errorbar_color=alex_orange)







        # T_on_dict = Dict(5 => 2.3364764812556524, 7 => 2.3147681253733903, 9 => 2.2833543431575998, 11 => 2.450317799353838) # Beta_stretch = 0.8*Beta(T=\infty)
        # epsilon_on_dict = Dict(5 => -0.21757, 7 => -0.2258468253968254, 11 => -0.2262537878787879, 9 => -0.2291601851851852)
    
        T_on_dict = Dict(5 => 1.9742874104754167, 7 => 1.9243708176391863, 9 => 2.018027471119676, 11 => 2.096487111125727) # Beta_stretch = 0.75*Beta(T=\infty)
        epsilon_on_dict = Dict(5 => -0.22990166666666664, 7 => -0.24041825396825398, 11 => -0.24010378787878783, 9 => -0.2408023148148148)

        # T_on_dict = Dict(5 => 1.728703038296018, 7 => 1.7614589684514228, 9 => 1.7819056391008394, 11 => 1.8886702059918599) # Beta_stretch = 0.7*Beta(T=\infty)
        # epsilon_on_dict = Dict(5 => -0.2422541666666667, 7 => -0.24992738095238093, 11 => -0.24862893939393932, 9 => -0.25057037037037044)
    
        # T_on_dict = Dict(5 => 1.6803133231596887, 7 => 1.6206595248088842, 9 => 1.6133644798131017, 11 => 1.6905658816913272) # Beta_stretch = beta_infinity - exp(-1)
        # epsilon_on_dict = Dict(5 => -0.24569250000000004, 7 => -0.25866904761904763, 11 => -0.26410803030303037, 9 => -0.2627643518518519)    


        # Add epsilon_on values to inset
        epsilon_on_values = [epsilon_on_dict[L] for L in Ls]
        scatter!(epsilon_parameters_graph, [1/L for L in Ls], epsilon_on_values, color=alex_red, label=L"\bar\epsilon^{\rm on}", markersize=4)
        fit = curve_fit(linear_model, [1/L for L in Ls], epsilon_on_values, [1.0, 0.0])
        m, c = fit.param
        confidence_intervals = confidence_interval(fit, 0.05)
        c_ci = confidence_intervals[2]
        c_error = c - c_ci[1]
        fit_1_L_values = LinRange(0.0, 0.2, 100)
        plot!(epsilon_parameters_graph, fit_1_L_values, m*fit_1_L_values .+ c, color=alex_red, label="", linestyle=:dash)
        scatter!(epsilon_parameters_graph, [0.0], [c], yerror=(c_error),color=alex_red, markerstrokecolor=alex_red, markersize=1, label="", errorbar_color=alex_red)



    end




    


    ### --- SAVE GRAPH ---

    savefig(graph, "results/relaxed_anneal_results/$(output_name).png")
    savefig(graph, "results/relaxed_anneal_results/$(output_name).pdf")
    display(graph)

    if model == "inherent_disorder"
        savefig(epsilon_parameters_graph, "results/relaxed_anneal_results/$(output_name)_epsilon_parameters.png")
        savefig(epsilon_parameters_graph, "results/relaxed_anneal_results/$(output_name)_epsilon_parameters.pdf")
        display(epsilon_parameters_graph)
    end

    ### --- PRINT ERRORS ---
    println("The following files do not exist:")
    for filename in filenames_that_do_not_exist
        println(filename)
    end

end














function temperature_to_energy_density_extractor(L::Int64, T::Float64; model::String="inherent_disorder", swap_move_probability::Float64=1.0)
 
    ## --- All L Figure ---
    models = [model]
    Ls = [L]
    swap_move_probabilities = [swap_move_probability]
    trials = 50
    N_T = 200


    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]
    custom_epsilon_star = []
    custom_epsilon_0 = []

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                if L==11 && swap_move_probability == 1.0 && model == "inherent_disorder"
                    N_T = 100
                else
                    N_T = 200
                end


                temperatures = zeros(N_T)
                running_total_average_normalised_energy_densities_by_temperature = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    # filename = "results/final_paper_results/relaxed_anneal_results/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    filename = "results/relaxed_anneal_results/data/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_average_normalised_energy_densities_by_temperature .+= data_matrix[:,3]
                        
                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_average_normalised_energy_densities_by_temperature/actual_number_of_trials

            end
        end
    end


    temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
    energy_density = results_dictionary[(model, L, swap_move_probability, "average_energy_densities")][argmin(abs.(temperatures.-T))]

 
    return energy_density
end




function energy_density_to_temperature_extractor(L::Int64, epsilon::Float64; model::String="inherent_disorder", swap_move_probability::Float64=1.0)
 
    ## --- All L Figure ---
    models = [model]
    Ls = [L]
    swap_move_probabilities = [swap_move_probability]
    trials = 50
    N_T = 200


    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                if L==11 && swap_move_probability == 1.0 && model == "inherent_disorder"
                    N_T = 100
                else
                    N_T = 200
                end


                temperatures = zeros(N_T)
                running_total_average_normalised_energy_densities_by_temperature = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    # filename = "results/final_paper_results/relaxed_anneal_results/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    filename = "results/relaxed_anneal_results/data/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_average_normalised_energy_densities_by_temperature .+= data_matrix[:,3]

                        
                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_average_normalised_energy_densities_by_temperature/actual_number_of_trials

            end
        end
    end

    energy_densities = results_dictionary[(model, L, swap_move_probability, "average_energy_densities")]
    temperature = results_dictionary[(model, L, swap_move_probability, "temperatures")][argmin(abs.(energy_densities.-epsilon))]



    return temperature
end

relaxed_anneal_figure()

all_L_relaxed_anneal_figure("inherent_disorder")

all_L_relaxed_anneal_figure("clean")