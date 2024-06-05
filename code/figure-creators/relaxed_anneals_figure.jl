using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using Images

include("../core/rubiks_cube.jl")

function relaxed_anneal_figure(output_name::String)

    ## --- All L Figure ---
    # models = ["clean", "inherent_disorder"]
    # Ls = [3, 5, 7, 9, 11]
    # swap_move_probabilities = [0.0, 1.0]
    # trials = 50
    # N_T = 100

    # --- L=11 Figure ---
    models = ["clean", "inherent_disorder"]
    Ls = [11]
    swap_move_probabilities = [0.0, 1.0]
    trials = 50

    ## --- Custom FIgure ---
    # models = ["custom"]
    # Ls = [11]
    # swap_move_probabilities = [0.0, 1.0]
    # trials = 50


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
    custom_epsilon_star = []
    custom_epsilon_0 = []

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                N_T = model == "clean" ? 200 : 100

                temperatures = zeros(N_T)
                running_total_average_energy_densities_by_temperature = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    filename = "results/final_paper_results/relaxed_anneal_results/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
                    try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_average_energy_densities_by_temperature .+= data_matrix[:,3]

                        if model=="custom"
                            println("$(model) $(L) $(swap_move_probability) $(minimum(data_matrix[:,3]))")
                            if swap_move_probability == 0.0
                                push!(custom_epsilon_star, minimum(data_matrix[:,3]))
                            elseif swap_move_probability == 1.0
                                push!(custom_epsilon_0, minimum(data_matrix[:,3]))
                            end

                        end
                        
                        actual_number_of_trials += 1
                    catch e
                        push!(filenames_that_do_not_exist, filename)
                    end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_average_energy_densities_by_temperature/actual_number_of_trials
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
        elseif model=="custom"
            colors = [alex_blue]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_energy_densities = results_dictionary[(model, L, swap_move_probability, "average_energy_densities")]

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

                plot!(graph,temperatures, average_energy_densities, label=label, color=colors[mod1(color_index,length(colors))], linestyle=linestyle)

            end
            color_index += 1
        end
    end

    ### --- ADD ANNOTATIONS TO GRAPH ---
    if "inherent_disorder" in models
        annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"\bar{\epsilon}^*", 12, alex_red, ))])
    elseif "custom" in models
        annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"{\epsilon}^*", 12, alex_blue, ))])
        # Add dot at (lowest T value with average_energy_densities values, custom average_energy_densities at T nearest to 0)
        custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
        scatter!([minimum(results_dictionary[("custom", 11, 0.0, "temperatures")])], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")]))]], color=alex_blue, label="", markersize=1.5)
    end

    if "inherent_disorder" in models
        T_star = 0.87
        annotate!(graph, [(T_star+0.35, ylims(graph)[2]-0.28, Plots.text(L"\bar{T}^*", 12, alex_red))])
    elseif "custom" in models
        T_star = 0.87
        annotate!(graph, [(T_star+0.4, ylims(graph)[2]-0.28, Plots.text(L"T^*", 12, alex_blue))])
        # Add dot at (T_star, custom average_energy_densities at T nearest to T_star)
        custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
        scatter!([T_star], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))]], color=alex_blue, label="", markersize=1.5)
    end


    if "clean" in models
        T_c = 0.90
        annotate!(graph, [(T_c+0.25, ylims(graph)[1]+0.05, Plots.text(L"T_c", 12, alex_alt_blue))])
    end


    if "custom" in models
        # Print epsilon^* for custom data if it exists (i.e. lowest value of epsilon for p_swap=0.0)
        custom_epsilon_star_avg = minimum(results_dictionary[("custom", 11, 0.0, "average_energy_densities")])
        println("Custom epsilon^* = $custom_epsilon_star_avg")

        # Print epsilon_0 for custom data if it exists (i.e. lowest value of epsilon for p_swap=1.0)
        custom_epsilon_0_avg = minimum(results_dictionary[("custom", 11, 1.0, "average_energy_densities")])
        println("Custom epsilon_0 = $custom_epsilon_0_avg")

        println("Mean epsilon^* = $(mean(custom_epsilon_star))")
        println("Standard Deviation epsilon^* = $(std(custom_epsilon_star))")

        println("Mean epsilon_0 = $(mean(custom_epsilon_0))")
        println("Standard Deviation epsilon_0 = $(std(custom_epsilon_0))")
    end

 

    ### --- SAVE GRAPH ---
    savefig(graph, "results/final_paper_results/$(output_name).png")
    savefig(graph, "results/final_paper_results/$(output_name).svg")
    display(graph)

    ### --- PRINT ERRORS ---
    println("The following files do not exist:")
    for filename in filenames_that_do_not_exist
        println(filename)
    end

end



   # ### --- PLOT ROTATIONS IMAGE ON GRAPH ---
    # img = load("results/final_paper_results/slice-rotations.png")

    # # Determine the desired width and height on the graph
    # # Here you set one dimension, and calculate the other to preserve the aspect ratio
    # desired_width = 0.55
    # aspect_ratio = size(img, 2) / size(img, 1) # width / height
    # desired_height = (desired_width / aspect_ratio)

    # # Determine the location on the graph where you want the image's bottom-left corner
    # x_location = 0.26
    # y_location = 0.21

    # # Calculate x and y ranges for the image placement
    # xrange = [x_location, x_location - desired_width]
    # yrange = [y_location, y_location - desired_height]

    # # Plot the image with the specified dimensions and location
    # # plot!(graph, xrange, yrange, reverse(img; dims=1), yflip=false, inset=bbox(x_location,y_location,desired_width,desired_height), subplot=2, aspect_ratio=:auto, axis=false, grid=false, framestyle=:none, legend=false, ticks=nothing, border=:none, plot_bgcolor=:transparent)

