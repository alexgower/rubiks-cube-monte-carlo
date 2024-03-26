using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

include("../core/rubiks_cube.jl")

function dynamic_autocorrelation_function_figure(simulation_name::String; type::String="configuration", display_temperatures::Vector{Float64}=empty([0.0]))


    ### --- SET UP DEFAULT PARAMETERS ---
    filename = "results/final_paper_results",simulation_name*"_$(type)_autocorrelation_averages_by_time.csv"
    header_line = readlines(joinpath(filename))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    dynamic_autocorrelation_function_data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)    

    ### -- CALCULATION ---
    sample_temperature = dynamic_autocorrelation_function_data_matrix[:,1]
    samples_in_average = dynamic_autocorrelation_function_data_matrix[:,2]
    dynamic_autocorrelation_functions_by_temperature = dynamic_autocorrelation_function_data_matrix[:,3:end]

    # For all repeated sample tempeartures, remove them from the sample_temperatures array 
    # And average (weighted by samples in each average) the corresponding dynamic autocorrelation functions
    unique_sample_temperatures = unique(sample_temperature)
    unique_dynamic_autocorrelation_functions = zeros(Float64, length(unique_sample_temperatures), size(dynamic_autocorrelation_functions_by_temperature,2))
    for (i, temperature) in pairs(unique_sample_temperatures)
        indices = findall(sample_temperature .== temperature)
        unique_dynamic_autocorrelation_functions[i,:] = sum(dynamic_autocorrelation_functions_by_temperature[indices,:].*samples_in_average[indices], dims=1) / sum(samples_in_average[indices])
    end


    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    Plots.default(dpi = 300)

    ### --- PLOTTING ---
    str = type == "configuration" ? "Configuration" : "Energy"
    symbol = type == "configuration" ? L"\mathcal{C}(t)" : L"\mathcal{C}_E(t)"
    graph = plot(title="", xlabel="Time [MC Steps]", ylabel="$(str) Autocorrelation Function, $(symbol)", legend=:topright)

    colors = [alex_red, alex_pink, alex_orange, alex_green, alex_blue]

    for (i, temperature) in pairs(unique_sample_temperatures)
        if isempty(display_temperatures) || temperature in unique_sample_temperatures
            plot!(graph, unique_dynamic_autocorrelation_functions[i,:], label="T = "*string(temperature), color=colors[mod1(i, length(colors))])
        end
    end

    ### --- SAVING GRAPHS ---
    savefig(graph, joinpath("results/final_paper_results",simulation_name*"_$(type)_autocorrelation_averages_by_time.png"))
    savefig(graph, joinpath("results/final_paper_results",simulation_name*"_$(type)_autocorrelation_averages_by_time.svg"))
    display(graph)




end