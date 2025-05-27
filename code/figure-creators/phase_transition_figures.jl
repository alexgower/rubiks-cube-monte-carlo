using Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using LsqFit

include("../core/rubiks_cube.jl")

gaussian(x, p) = p[1] * exp.(-(x .- p[2]).^2 ./ (2 * p[3]^2)) .+ p[4]

function phase_transition_figures()

    # --- All L Figure ---
    # models = ["clean", "inherent_disorder"]
    models = ["clean","inherent_disorder"]
    Ls = [11,9,7,5,3]
    swap_move_probabilities = [1.0]
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
                running_total_specific_heat_capacities = zeros(N_T)
                running_total_specific_heat_capacities_squared = zeros(N_T)  # For error calculation
                running_total_binder_cumulants = zeros(N_T)
                running_total_binder_cumulants_squared = zeros(N_T)  # For error calculation

                actual_number_of_trials=0
                for trial in 1:trials
                    # filename = "results/final_paper_results/relaxed_anneal_results/" * model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    filename = "results/relaxed_anneal_results/data/" * model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"


                    try
                        println("Reading in data from: ", filename)
                        println("Using N_T = $(N_T)")
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        specific_heat_capacities = data_matrix[:,5] ./ (6*L^2)
                        running_total_specific_heat_capacities .+= specific_heat_capacities
                        running_total_specific_heat_capacities_squared .+= specific_heat_capacities.^2

                        # 7 = <M^2>, 8=<M^4>
                        binder_cumulants = 1 .- (data_matrix[:,8] ./ (3 * data_matrix[:,7].^2))
                        running_total_binder_cumulants .+= binder_cumulants
                        running_total_binder_cumulants_squared .+= binder_cumulants.^2

                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                end

                println("Actual number of trials for $(model) L = $L: ", actual_number_of_trials)

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures


                # Calculate mean and error for specific heat
                mean_specific_heat = running_total_specific_heat_capacities/actual_number_of_trials
                mean_specific_heat_squared = running_total_specific_heat_capacities_squared/actual_number_of_trials
                specific_heat_variance = mean_specific_heat_squared .- mean_specific_heat.^2
                specific_heat_error = sqrt.(max.(0, specific_heat_variance)/actual_number_of_trials)

                results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")] = mean_specific_heat
                results_dictionary[(model, L, swap_move_probability, "specific_heat_errors")] = specific_heat_error

                # Calculate mean and error for Binder cumulants
                mean_binder = running_total_binder_cumulants/actual_number_of_trials
                mean_binder_squared = running_total_binder_cumulants_squared/actual_number_of_trials
                binder_variance = mean_binder_squared .- mean_binder.^2
                binder_error = sqrt.(max.(0, binder_variance)/actual_number_of_trials)

                results_dictionary[(model, L, swap_move_probability, "average_binder_cumulants")] = mean_binder
                results_dictionary[(model, L, swap_move_probability, "binder_cumulant_errors")] = binder_error
            end
        end
    end



    ### --- PLOT SPECIFIC HEAT CAPACITIES DATA ---
    for model in models
        # Set different ylims based on model
        ylimit = model == "clean" ? 8.5 : 2.0
        graph = plot(title="", xlabel="Temperature, "*L"T", legend=:bottomright, yaxis="Specific Heat Capacity, "*L"\langle c \rangle = \frac{\langle  E^2 \rangle - \langle\! E \rangle^2}{6 L^2 T^2}", margin=3mm, ylabelfontsize=10, xlabelfontsize=10, ylims=(0,ylimit))

        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        elseif model=="custom"
            colors = [alex_blue]
        end

        for (i, L) in enumerate(Ls)
            if L==3
                continue
            end

            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_specific_heat_capacities = results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")]
                specific_heat_errors = results_dictionary[(model, L, swap_move_probability, "specific_heat_errors")]

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

                plot!(graph, temperatures, average_specific_heat_capacities, 
                ribbon=specific_heat_errors,
                fillalpha=0.2,
                label=label, 
                color=colors[mod1(i,length(colors))], 
                linestyle=linestyle)

                # Also plot a dashed line at the height of the maximum
                # max_specific_heat_capacity = maximum(average_specific_heat_capacities)
                # hline!([max_specific_heat_capacity], color=colors[mod1(color_index,length(colors))], linestyle=:dash, label="")
            end
        end


        # # Plot an inset with the maximum specific heat capacity per L against L
        # max_specific_heat_capacities = [maximum(results_dictionary[(model, L, 1.0, "average_specific_heat_capacities")]) for L in Ls]
        # # plot!(graph, Ls, max_specific_heat_capacities, label="", inset=bbox(0.4,0.15,0.3,0.4), subplot=2, xlabel=L"L", ylabel=L"\langle C \rangle_{\rm{max}}", yguidefontsize=12,xguidefontsize=12, color=colors[mod1(color_index,length(colors))])
        # scatter!(graph, Ls.^2, max_specific_heat_capacities, label="", inset=bbox(0.4,0.15,0.3,0.4), subplot=2, xlabel=L"L^2", ylabel=L"\langle C \rangle_{\rm{max}}", yguidefontsize=12,xguidefontsize=12, color=colors[mod1(color_index,length(colors))])


        max_specific_heat_capacities = Float64[]
        max_specific_heat_capacities_errors = Float64[]  # New array for errors
        max_temperatures = Float64[]
        gaussian_widths = Float64[]  # New array for standard deviations
        gaussian_widths_errors = Float64[]  # New array for standard deviation errors

        for L in Ls
            temperatures = results_dictionary[(model, L, 1.0, "temperatures")]
            average_specific_heat_capacities = results_dictionary[(model, L, 1.0, "average_specific_heat_capacities")]

            println("Average specific heat capacities for L = $L: ", average_specific_heat_capacities)
            
            max_c, max_t, max_c_err, width, width_err = fit_gaussian_peak(temperatures, average_specific_heat_capacities, model)
            push!(max_specific_heat_capacities, max_c)
            push!(max_specific_heat_capacities_errors, max_c_err)  # Store the error
            push!(max_temperatures, max_t)
            push!(gaussian_widths, width)  # Store the width
            push!(gaussian_widths_errors, width_err)  # Store the width error

            println("L = $L, max_c = $max_c ± $max_c_err, max_T = $max_t, width = $width ± $width_err")
        end

        # Create first inset plot (peak heights)
        plot!(graph, [], [], inset=bbox(0.32,0.05,0.3,0.3), subplot=2, 
              xlabel=L"L^2", ylabel=L"\langle c \rangle_{\rm{max}}", 
              yguidefontsize=12, xguidefontsize=12, legend=:none)
        
        # Create second inset plot (widths)
        yticks = model == "clean" ? [0.01,0.015,0.02,0.025,0.03] : [0.07,0.08,0.09,0.10,0.11]
        plot!(graph, [], [], inset=bbox(0.32,0.45,0.3,0.3), subplot=3, 
              xlabel=L"L^2", ylabel=L"\sigma", 
              yguidefontsize=12, xguidefontsize=12, legend=:none, yticks=yticks)
        
        # Add points with matching colors and error bars to both insets
        for (i, L) in enumerate(Ls)
            if L==3
                continue
            end
            # Plot in first inset (peak heights)
            scatter!(graph, [L^2], [max_specific_heat_capacities[i]], 
                    yerror=max_specific_heat_capacities_errors[i],
                    subplot=2, 
                    color=colors[mod1(i,length(colors))], 
                    label="")
            
            # Plot in second inset (widths) with error bars
            scatter!(graph, [L^2], [gaussian_widths[i]], 
                    yerror=gaussian_widths_errors[i],
                    subplot=3, 
                    color=colors[mod1(i,length(colors))], 
                    label="")
        end

        # Add linear fit to the insets only for clean model
        if model == "clean"
            # Filter out L=3 data point
            valid_indices = findall(x -> x != 3, Ls)
            x_data = Float64[Ls[i]^2 for i in valid_indices]
            
            # Fit for peak heights
            y_data = max_specific_heat_capacities[valid_indices]
            y_errors = max_specific_heat_capacities_errors[valid_indices]
            
            # Weighted linear fit using error bars
            weights = 1.0 ./ (y_errors.^2)
            X = hcat(ones(length(x_data)), x_data)
            W = Diagonal(weights)
            β = (X' * W * X) \ (X' * W * y_data)
            
            # Calculate fit line points
            x_fit = range(minimum(x_data), maximum(x_data), length=100)
            y_fit = β[1] .+ β[2] .* x_fit
            
            # Plot the fit line for peak heights
            plot!(graph, x_fit, y_fit, 
                  subplot=2,
                  color=:black,
                  linestyle=:dash,
                  label="",
                  linewidth=1)
            
            # Print fit parameters
            println("\nLinear fit parameters for clean model (peak heights):")
            println("Intercept: $(round(β[1], digits=4))")
            println("Slope: $(round(β[2], digits=4))")
            
            # Fit for widths
            y_data = gaussian_widths[valid_indices]
            
            # Simple linear fit for widths (no error bars)
            X = hcat(ones(length(x_data)), x_data)
            β = X \ y_data
            
            # Calculate fit line points
            y_fit = β[1] .+ β[2] .* x_fit
            
            # Plot the fit line for widths
            plot!(graph, x_fit, y_fit, 
                  subplot=3,
                  color=:black,
                  linestyle=:dash,
                  label="",
                  linewidth=1)
            
            println("\nLinear fit parameters for clean model (widths):")
            println("Intercept: $(round(β[1], digits=4))")
            println("Slope: $(round(β[2], digits=4))")
        end

        ### --- SAVE GRAPH --- 
        # savefig(graph, "results/final_paper_results/specific_heat_capacity_L_scaling_$(model).png")
        # savefig(graph, "results/final_paper_results/specific_heat_capacity_L_scaling_$(model).pdf")
        savefig(graph, "results/relaxed_anneal_results/specific_heat_capacity_L_scaling_$(model)_inset.png")
        savefig(graph, "results/relaxed_anneal_results/specific_heat_capacity_L_scaling_$(model)_inset.pdf")
        display(graph)
    end










    ### --- PLOT BINDER CUMULANTS DATA ---
    for model in models
        graph = plot(title="", xlabel="Temperature, "*L"T", legend=:topright, yaxis="Binder Cumulant, "*L"U_L = 1 - \frac{\langle M^4 \rangle_L}{3\langle M^2 \rangle_L^2}", margin=3mm, ylabelfontsize=10, xlabelfontsize=10)


        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        elseif model=="custom"
            colors = [alex_blue]
        end

        for (i, L) in enumerate(Ls)
            if L==3
                continue
            end
            
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_binder_cumulants = results_dictionary[(model, L, swap_move_probability, "average_binder_cumulants")]
                binder_cumulant_errors = results_dictionary[(model, L, swap_move_probability, "binder_cumulant_errors")]

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

                # Plot with error bars
                plot!(graph, temperatures, average_binder_cumulants, 
                      ribbon=binder_cumulant_errors,
                      fillalpha=0.2,
                      label=label, 
                      color=colors[mod1(i,length(colors))], 
                      linestyle=linestyle)
            end 
        end

        ### --- SAVE GRAPH ---
        savefig(graph, "results/relaxed_anneal_results/binder_cumulant_L_scaling_$(model).png")
        savefig(graph, "results/relaxed_anneal_results/binder_cumulant_L_scaling_$(model).pdf")
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

















    # TODO remove
    ### --- Analyze peak specific heat for all L ---
    println("\nAnalyzing peak specific heat capacity for all L:")
    for model in models
        for L in Ls
            if L==3  # Skip L=3 as we do in main plots
                continue
            end
            
            swap_move_probability = 1.0
            
            # First find the temperature index where average is maximum
            temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
            average_specific_heat = results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")]
            peak_index = argmax(average_specific_heat)
            peak_temperature = temperatures[peak_index]
            
            println("\nModel: $(model), L = $(L)")
            println("Peak temperature: $(peak_temperature)")
            println("Average specific heat at peak: $(average_specific_heat[peak_index])")
            println("\nIndividual trial values at T=$(peak_temperature):")
            
            # Create plot for all trials for this L
            trial_plot = plot(title="L=$(L) Specific Heat Capacity - All Trials ($(model))", 
                            xlabel="Temperature, "*L"T", 
                            ylabel="Specific Heat Capacity, "*L"\langle c \rangle",
                            legend=:none)
            
            # Now collect all trial values at this temperature and plot full curves
            peak_values = Float64[]
            for trial in 1:trials
                filename = "results/relaxed_anneal_results/data/" * model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                try
                    data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                    specific_heat = data_matrix[:,5] ./ (6*L^2)
                    push!(peak_values, specific_heat[peak_index])
                    println("Trial $(trial): $(specific_heat[peak_index])")
                    
                    # Plot individual trial curve with low alpha
                    plot!(trial_plot, temperatures, specific_heat, 
                          color=:grey, 
                          alpha=0.3, 
                          linewidth=1)
                catch e
                    continue
                end
            end
            
            # Plot average curve on top with error ribbon
            specific_heat_errors = results_dictionary[(model, L, swap_move_probability, "specific_heat_errors")]
            plot!(trial_plot, temperatures, average_specific_heat,
                  ribbon=specific_heat_errors,
                  color=model == "clean" ? alex_alt_blue : alex_red,
                  linewidth=2,
                  label="Average")
            
            println("\nStatistics at peak:")
            println("Mean: $(mean(peak_values))")
            println("Standard deviation: $(std(peak_values))")
            println("Standard error: $(std(peak_values)/sqrt(length(peak_values)))")
            
            # Save the trial plot
            savefig(trial_plot, "results/relaxed_anneal_results/rough/L$(L)_all_trials_$(model).png")
            savefig(trial_plot, "results/relaxed_anneal_results/rough/L$(L)_all_trials_$(model).pdf")
            display(trial_plot)
        end
    end

    ### --- Analyze Gaussian peak fitting ---
    println("\nAnalyzing Gaussian peak fitting for all L:")
    for model in models
        for L in Ls
            if L==3  # Skip L=3 as we do in main plots
                continue
            end
            
            swap_move_probability = 1.0
            temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
            average_specific_heat = results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")]
            
            # Create plot for Gaussian fit
            fit_plot = plot(title="L=$(L) Specific Heat Capacity - Gaussian Fit ($(model))", 
                          xlabel="Temperature, "*L"T", 
                          ylabel="Specific Heat Capacity, "*L"\langle c \rangle",
                          legend=:topright)
            
            # Plot the data
            plot!(fit_plot, temperatures, average_specific_heat,
                  color=model == "clean" ? alex_alt_blue : alex_red,
                  linewidth=2,
                  label="Data")
            
            # Get peak values using existing function
            peak_height, peak_temp, peak_height_error, width = fit_gaussian_peak(temperatures, average_specific_heat, model)
            
            # Generate fitted curve for visualization
            x_fit = LinRange(minimum(temperatures), maximum(temperatures), 1000)
            x_range = temperatures[end] - temperatures[1]
            fit_fraction = model == "clean" ? 0.005 : 0.02
            fit_range = (temperatures .>= peak_temp - fit_fraction*x_range) .& (temperatures .<= peak_temp + fit_fraction*x_range)
            
            # Reuse the gaussian function from fit_gaussian_peak for visualization
            p0 = [peak_height - minimum(average_specific_heat), peak_temp, fit_fraction*x_range, minimum(average_specific_heat)]
            fit = curve_fit(gaussian, temperatures[fit_range], average_specific_heat[fit_range], p0)
            y_fit = gaussian(x_fit, fit.param)
            
            # Calculate actual peak height from fit parameters
            fitted_peak_height = fit.param[1] + fit.param[4]  # amplitude + offset
            
            # Plot the fit
            plot!(fit_plot, x_fit, y_fit,
                  color=:pink,
                #   linestyle=,
                  linewidth=2,
                  label="Gaussian Fit")
            
            # Add vertical line at peak
            vline!(fit_plot, [peak_temp],
                  color=:gray,
                  linestyle=:dot,
                  label="Peak T = $(round(peak_temp, digits=3))")
            
            # Add horizontal line at fitted peak height
            hline!(fit_plot, [fitted_peak_height],
                  color=:gray,
                  linestyle=:dash,
                  label="Peak Height = $(round(fitted_peak_height, digits=3))")
            
            # Print fit parameters
            println("\nModel: $(model), L = $(L)")
            println("Gaussian fit parameters:")
            println("Peak position = $(round(peak_temp, digits=4))")
            println("Peak height = $(round(fitted_peak_height, digits=4))")
            println("Standard deviation = $(round(width, digits=4))")
            println("Baseline offset = $(round(fit.param[4], digits=4))")
            
            # Save the fit plot
            savefig(fit_plot, "results/relaxed_anneal_results/rough/L$(L)_gaussian_fit_$(model).png")
            savefig(fit_plot, "results/relaxed_anneal_results/rough/L$(L)_gaussian_fit_$(model).pdf")
            display(fit_plot)
        end
    end
end



function fit_gaussian_peak(x, y, model)
    gaussian(x, p) = p[1] * exp.(-(x .- p[2]).^2 ./ (2 * p[3]^2)) .+ p[4]

    fit_fraction = model == "clean" ? 0.005 : 0.02

    # Find the initial guess for the peak
    max_index = argmax(y)
    max_y = y[max_index]
    max_x = x[max_index]
    x_range = abs(x[end] - x[1])
    
    # Initial parameters: [amplitude, mean, std_dev, offset]
    p0 = [max_y - minimum(y), max_x, x_range*fit_fraction, minimum(y)]
    
    # Fit only around the peak
    fit_range = (x .>= max_x - fit_fraction*x_range) .& (x .<= max_x + fit_fraction*x_range)
    x_fit = x[fit_range]
    y_fit = y[fit_range]

    println("Max y: $(max_y)")
    println("Max x: $(max_x)")
    println("x range: $(x_range)")
    println("Fit fraction: $(fit_fraction)")
    println("Length of fit range: $(length(x[fit_range]))")

    # Perform the fit
    fit = curve_fit(gaussian, x[fit_range], y[fit_range], p0)
    
    # Extract the fitted parameters
    amplitude, mean, std_dev, offset = fit.param
    
    # Get parameter errors from the covariance matrix
    param_errors = stderror(fit)
    # Error propagation for peak height (amplitude + offset)
    peak_height_error = sqrt(param_errors[1]^2 + param_errors[4]^2)
    
    # Return peak height, position, error in peak height, and width
    return amplitude + offset, mean, peak_height_error, std_dev, param_errors[3]
end

phase_transition_figures()