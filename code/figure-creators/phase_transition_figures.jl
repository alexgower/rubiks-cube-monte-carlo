using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using LsqFit

include("../core/rubiks_cube.jl")

function phase_transition_figures()

    # --- All L Figure ---
    # models = ["clean", "inherent_disorder"]
    models = ["clean","inherent_disorder"]
    Ls = [11, 9,7,5,3]
    swap_move_probabilities = [1.0]
    trials = 50
    N_T = 100


    ### --- COLOURS ---
    Plots.default(dpi = 300)

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

                temperatures = zeros(N_T)
                running_total_specific_heat_capacities = zeros(N_T)
                running_total_binder_cumulants = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    filename = "results/final_paper_results/relaxed_anneal_results/" * model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    if model == "clean" && L==11
                        filename = filename * "_old"
                    end

                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_specific_heat_capacities .+= data_matrix[:,5] 
                        # 7 = <M^2>, 8=<M^4>
                        running_total_binder_cumulants .+= 1 .- (data_matrix[:,8] ./ (3 * data_matrix[:,7].^2))

                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")] = running_total_specific_heat_capacities/actual_number_of_trials
                results_dictionary[(model, L, swap_move_probability, "average_binder_cumulants")] = running_total_binder_cumulants/actual_number_of_trials
            end
        end
    end



    ### --- PLOT SPECIFIC HEAT CAPACITIES DATA ---
    color_index = 1
    for model in models
        graph = plot(title="", xlabel="Temperature, "*L"T", legend=:bottomright, yaxis="Heat Capacity, "*L"\langle C \rangle = \frac{\langle  E^2 \rangle - \langle\! E \rangle^2}{T^2}", margin=3mm, ylabelfontsize=10, xlabelfontsize=10)


        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        elseif model=="custom"
            colors = [alex_blue]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_specific_heat_capacities = results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")]

                linestyle = swap_move_probability == 0.0 ? :dash : :solid

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

                plot!(graph,temperatures, average_specific_heat_capacities, label=label, color=colors[mod1(color_index,length(colors))], linestyle=linestyle)

                # Also plot a dashed line at the height of the maximum
                # max_specific_heat_capacity = maximum(average_specific_heat_capacities)
                # hline!([max_specific_heat_capacity], color=colors[mod1(color_index,length(colors))], linestyle=:dash, label="")
            end
            color_index += 1
        end


        # # Plot an inset with the maximum specific heat capacity per L against L
        # max_specific_heat_capacities = [maximum(results_dictionary[(model, L, 1.0, "average_specific_heat_capacities")]) for L in Ls]
        # # plot!(graph, Ls, max_specific_heat_capacities, label="", inset=bbox(0.4,0.15,0.3,0.4), subplot=2, xlabel=L"L", ylabel=L"\langle C \rangle_{\rm{max}}", yguidefontsize=12,xguidefontsize=12, color=colors[mod1(color_index,length(colors))])
        # scatter!(graph, Ls.^2, max_specific_heat_capacities, label="", inset=bbox(0.4,0.15,0.3,0.4), subplot=2, xlabel=L"L^2", ylabel=L"\langle C \rangle_{\rm{max}}", yguidefontsize=12,xguidefontsize=12, color=colors[mod1(color_index,length(colors))])


        max_specific_heat_capacities = Float64[]
        max_temperatures = Float64[]
        for L in Ls
            temperatures = results_dictionary[(model, L, 1.0, "temperatures")]
            average_specific_heat_capacities = results_dictionary[(model, L, 1.0, "average_specific_heat_capacities")]
            
            max_c, max_t = fit_gaussian_peak(temperatures, average_specific_heat_capacities)
            push!(max_specific_heat_capacities, max_c)
            push!(max_temperatures, max_t)

            println("L = $L, max_C = $max_c, max_T = $max_t")
        end
        scatter!(graph, Ls.^2, max_specific_heat_capacities, label="", inset=bbox(0.4,0.15,0.3,0.4), subplot=2, xlabel=L"L^2", ylabel=L"\langle C \rangle_{\rm{max}}", yguidefontsize=12,xguidefontsize=12, color=colors[mod1(color_index,length(colors))])


        ### --- SAVE GRAPH --- 
        savefig(graph, "results/final_paper_results/specific_heat_capacity_L_scaling_$(model).png")
        savefig(graph, "results/final_paper_results/specific_heat_capacity_L_scaling_$(model).pdf")
        display(graph)
    end










    ### --- PLOT BINDER CUMULANTS DATA ---
    color_index = 1
    for model in models
        graph = plot(title="", xlabel="Temperature, "*L"T", legend=:topright, yaxis="Binder Cumulant, "*L"U_L = 1 - \frac{\langle M^4 \rangle_L}{3\langle M^2 \rangle_L^2}", margin=3mm, ylabelfontsize=10, xlabelfontsize=10)


        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        elseif model=="custom"
            colors = [alex_blue]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_binder_cumulants = results_dictionary[(model, L, swap_move_probability, "average_binder_cumulants")]

                linestyle = swap_move_probability == 0.0 ? :dash : :solid

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

                plot!(graph,temperatures, average_binder_cumulants, label=label, color=colors[mod1(color_index,length(colors))], linestyle=linestyle)

            end
            color_index += 1
        end

        ### --- SAVE GRAPH ---
        savefig(graph, "results/final_paper_results/binder_cumulant_L_scaling_$(model).png")
        savefig(graph, "results/final_paper_results/binder_cumulant_L_scaling_$(model).pdf")
        display(graph)
    end

    ### --- ADD ANNOTATIONS TO GRAPH ---
    # if "inherent_disorder" in models
    #     annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"\bar{\epsilon}^*", 12, alex_red, ))])
    #     # Add dot at (lowest T value with average_energy_densities values, inherent_disorder average_energy_densities at T nearest to 0)
    #     # inherent_disorder_average_energy_densities = results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")]
    #     # scatter!([minimum(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")])], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")]))]], color=alex_red, label="", markersize=1.5)
    # elseif "custom" in models
    #     annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"{\epsilon}^*", 12, alex_blue, ))])
    #     # Add dot at (lowest T value with average_energy_densities values, custom average_energy_densities at T nearest to 0)
    #     custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
    #     scatter!([minimum(results_dictionary[("custom", 11, 0.0, "temperatures")])], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")]))]], color=alex_blue, label="", markersize=1.5)
    # end

    # if "inherent_disorder" in models
    #     T_star = 0.8
    #     annotate!(graph, [(T_star+0.2, ylims(graph)[2]-0.2, Plots.text(L"\bar{T}^*", 12, alex_red))])
    #     # Add dot at (T_star, inherent_disorder average_energy_densities at T nearest to T_star)
    #     # inherent_disorder_average_energy_densities = results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")]
    #     # scatter!([T_star], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_star))]], color=alex_red, label="", markersize=1.5)
    # elseif "custom" in models
    #     T_star = 0.87
    #     annotate!(graph, [(T_star+0.4, ylims(graph)[2]-0.28, Plots.text(L"T^*", 12, alex_blue))])
    #     # Add dot at (T_star, custom average_energy_densities at T nearest to T_star)
    #     custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
    #     scatter!([T_star], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))]], color=alex_blue, label="", markersize=1.5)
    # end


    ### --- PRINT ERRORS ---
    println("The following files do not exist:")
    for filename in filenames_that_do_not_exist
        println(filename)
    end

end



function fit_gaussian_peak(x, y)
    gaussian(x, p) = p[1] * exp.(-(x .- p[2]).^2 ./ (2 * p[3]^2)) .+ p[4]

    # Find the initial guess for the peak
    max_index = argmax(y)
    max_y = y[max_index]
    max_x = x[max_index]
    
    # Initial parameters: [amplitude, mean, std_dev, offset]
    p0 = [max_y - minimum(y), max_x, (x[end] - x[1])/10, minimum(y)]
    
    # Fit only around the peak (e.g., ±20% of the x range)
    x_range = x[end] - x[1]
    fit_range = (x .>= max_x - 0.2*x_range) .& (x .<= max_x + 0.2*x_range)
    
    # Perform the fit
    fit = curve_fit(gaussian, x[fit_range], y[fit_range], p0)
    
    # Extract the fitted parameters
    amplitude, mean, std_dev, offset = fit.param
    
    # The maximum of the Gaussian is at the mean
    return amplitude + offset, mean
end