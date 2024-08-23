using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using Images

include("../core/rubiks_cube.jl")

function relaxed_anneal_figure(output_name::String; extraction::Bool=false)

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

    ## --- Custom Figure ---
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

    if "clean" in models 
        # Print average epsilon at highest temperature value 
        # for L=11, swap_move_probability=1.0 and L=11, swap_move_probability=0.0

        println("Clean epsilon at highest temperature value")
        println("L=11, swap_move_probability=1.0: $(results_dictionary[("clean", 11, 1.0, "average_energy_densities")][1])")
        println("L=11, swap_move_probability=0.0: $(results_dictionary[("clean", 11, 0.0, "average_energy_densities")][1])")
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
        # annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"\bar{\epsilon}^*", 12, alex_red, ))])

        T_star = 0.79
        T_on = 2.25

        println("T* Inherent Disorder = $T_star")
        println("T_on Inherent Disorder = $T_on")
        println("epsilon at T* Inherent Disorder = $(results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_star))])")
        annotate!(graph, [(0.45, ylims(graph)[2]-0.2, Plots.text(L"(\bar{T}^{\!\!\!*}\!\!\!\!,\bar{\epsilon}^{\!\!\!*}\!\!)", 12, alex_red))])
        annotate!(graph, [(T_on + 0.2, ylims(graph)[2]-0.06, Plots.text(L"(\bar{T}^{on}\!\!\!\!,\bar{\epsilon}^{on})", 12, alex_red))])


        # Add dot at (T_star, inherent_disorder average_energy_densities at T nearest to T_star)
        inherent_disorder_average_energy_densities = results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")]
        scatter!([T_star], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_star))]], color=alex_red, label="", markersize=1.5)
                

        # Get \epsilon_on as the average energy density  at T_on
        inherent_disorder_epsilon_on = inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")].-T_on))]
        println("Inherent Disorder epsilon_on = $inherent_disorder_epsilon_on")        
        # Add dot at (T_on, inherent_disorder average_energy_densities at T nearest to T_on)
        scatter!([T_on], [inherent_disorder_epsilon_on], color=alex_red, label="", markersize=1.5)

        if extraction

            # Add dot at (lowest T value with average_energy_densities values, inherent_disorder average_energy_densities at T nearest to 0)
            # scatter!([minimum(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")])], [inherent_disorder_average_energy_densities[argmin(abs.(results_dictionary[("inherent_disorder", 11, 0.0, "temperatures")]))]], color=alex_red, label="", markersize=1.5)

            inherent_disorder_epsilon_star_avg = minimum(results_dictionary[("inherent_disorder", 11, 0.0, "average_energy_densities")])
            inherent_disorder_epsilon_0_avg = minimum(results_dictionary[("inherent_disorder", 11, 1.0, "average_energy_densities")])

            println("Inherent Disorder epsilon^* = $inherent_disorder_epsilon_star_avg")
            println("Inherent Disorder epsilon_0 = $inherent_disorder_epsilon_0_avg")



        end
    end

    if "clean" in models
        if extraction
            annotate!(graph, [(0.3, ylims(graph)[2]-0.35, Plots.text(L"{\epsilon}^*", 12, alex_alt_blue, ))])
        end

        T_c = 0.9
        println("T_c = $T_c")
        annotate!(graph, [(T_c+0.25, ylims(graph)[1]+0.05, Plots.text(L"T_c", 12, alex_alt_blue))])

        T_star = 0.92
        println("T* Clean = $T_star")
        println("epsilon at T* Clean = $(results_dictionary[("clean", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")].-T_star))])")
        if extraction
            annotate!(graph, [(T_star+0.6, ylims(graph)[2]-0.3, Plots.text(L"T^*", 12, alex_alt_blue))])
        end

        if extraction
            # Add dot at (T_star, clean_average_energy_densities at T nearest to T_star)
            clean_average_energy_densities = results_dictionary[("clean", 11, 0.0, "average_energy_densities")]
            scatter!([T_star], [clean_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")].-T_star))]], color=alex_alt_blue, label="", markersize=1.5)
        
            # Add dot at (lowest T value with average_energy_densities values, clean_average_energy_densities at T nearest to 0)
            scatter!([minimum(results_dictionary[("clean", 11, 0.0, "temperatures")])], [clean_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 0.0, "temperatures")]))]], color=alex_alt_blue, label="", markersize=1.5)

            # Add dot at (T_c, clean-average_energy_densities at T nearest to T_c)
            clean_swap_average_energy_densities = results_dictionary[("clean", 11, 1.0, "average_energy_densities")]
            scatter!([T_c], [clean_swap_average_energy_densities[argmin(abs.(results_dictionary[("clean", 11, 1.0, "temperatures")].-T_c))]], color=alex_alt_blue, label="", markersize=1.5)

            clean_epsilon_star_avg = minimum(results_dictionary[("clean", 11, 0.0, "average_energy_densities")])
            clean_epsilon_0_avg = minimum(results_dictionary[("clean", 11, 1.0, "average_energy_densities")])

            println("Clean epsilon^* = $clean_epsilon_star_avg")
            println("Clean epsilon_0 = $clean_epsilon_0_avg")
        end
    end

    if "custom" in models
        annotate!(graph, [(0.35, ylims(graph)[2]-0.24, Plots.text(L"{\epsilon}^*", 12, alex_blue, ))])

        T_star = 0.86
        println("T* Custom = $T_star")
        println("epsilon at T* Custom = $(results_dictionary[("custom", 11, 0.0, "average_energy_densities")][argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))])")
        annotate!(graph, [(T_star+0.4, ylims(graph)[2]-0.28, Plots.text(L"T^*", 12, alex_blue))])
        if extraction
            # Add dot at (T_star, custom average_energy_densities at T nearest to T_star)
            custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
            scatter!([T_star], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")].-T_star))]], color=alex_blue, label="", markersize=1.5)

            # Add dot at (lowest T value with average_energy_densities values, custom average_energy_densities at T nearest to 0)
            custom_average_energy_densities = results_dictionary[("custom", 11, 0.0, "average_energy_densities")]
            scatter!([minimum(results_dictionary[("custom", 11, 0.0, "temperatures")])], [custom_average_energy_densities[argmin(abs.(results_dictionary[("custom", 11, 0.0, "temperatures")]))]], color=alex_blue, label="", markersize=1.5)
        
                       
            custom_epsilon_star_avg = minimum(results_dictionary[("custom", 11, 0.0, "average_energy_densities")])
            custom_epsilon_0_avg = minimum(results_dictionary[("custom", 11, 1.0, "average_energy_densities")])

            println("Custom epsilon^* = $custom_epsilon_star_avg")
            println("Custom epsilon_0 = $custom_epsilon_0_avg")
        end
    end



 

    ### --- SAVE GRAPH ---
    savefig(graph, "results/final_paper_results/$(output_name).png")
    savefig(graph, "results/final_paper_results/$(output_name).pdf")
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

