using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using Images

include("../core/rubiks_cube.jl")

function relaxed_anneal_figures()

    ### --- DATA TO INCLUDE ---
    # models = ["clean", "inherent_disorder", "custom"]
    # Ls = [3, 5, 7, 9, 11]
    # swap_move_probabilities = [0.0, 1.0]
    # trials = 50
    # N_T = 100

    models = ["inherent_disorder"]
    Ls = [11]
    swap_move_probabilities = [0.0, 1.0]
    trials = 3
    N_T = 99

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    Plots.default(dpi = 300)

    ### --- READ IN DATA ---
    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                temperatures = zeros(N_T)
                running_total_average_energy_densities_by_temperature = zeros(N_T)

                for trial in 1:trials
                    filename = "results/final_paper_results/relaxed_anneal_results/", model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"

                    data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                    
                    temperatures .= data_matrix[:,1]
                    running_total_average_energy_densities_by_temperature .+= data_matrix[:,3]
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energy_densities")] = running_total_average_energy_densities_by_temperature/trials
            end
        end
    end


    ### --- PLOT DATA ---
    graph = plot(title="", xlabel="Temperature, "*L"T", legend=:bottomright, yaxis="Average Energy Density, "*L"\langle\! E/|\!\!E_s|\!\rangle", ylims=(-1.0,-0.1))

    colors = [alex_red, alex_pink, alex_orange, alex_green, alex_blue]

    color_index = 1
    for model in models
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
                end

                plot!(graph,temperatures, average_energy_densities, label=label, color=colors[color_index], linestyle=linestyle)

                color_index += 1
            end
        end
    end

    ### --- PLOT IMAGE ON GRAPH ---
    img = load("results/final_paper_results/slice-rotations.png")

    # Determine the desired width and height on the graph
    # Here you set one dimension, and calculate the other to preserve the aspect ratio
    desired_width = 0.4
    aspect_ratio = size(img, 2) / size(img, 1) # width / height
    desired_height = (desired_width / aspect_ratio)

    # Determine the location on the graph where you want the image's bottom-left corner
    x_location = 0.25
    y_location = 0.4

    # Calculate x and y ranges for the image placement
    xrange = [x_location, x_location - desired_width]
    yrange = [y_location, y_location - desired_height]

    # Plot the image with the specified dimensions and location
    plot!(graph, xrange, yrange, reverse(img; dims=1), yflip=false, inset=bbox(x_location,y_location,desired_width,desired_height), subplot=2, aspect_ratio=:auto, axis=false, grid=false, framestyle=:none, legend=false, ticks=nothing, border=:none, plot_bgcolor=:transparent)



    ### --- SAVE GRAPH ---
    savefig(graph, "results/final_paper_results/relaxed_anneals_figure.png")
    savefig(graph, "results/final_paper_results/relaxed_anneals_figure.svg")
    display(graph)

end

# TODO
# - make error bars graphs with average energies squared
# - make specific heat capacity graphs