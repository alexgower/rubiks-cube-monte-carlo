using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

include("../core/rubiks_cube.jl")

function stretched_exponential_parameters_figure(simulation_name::String)


    ### --- SET UP DEFAULT PARAMETERS ---
    filename = "results/final_paper_results",simulation_name*"_configuration_taus.csv"
    header_line = readlines(joinpath(filename))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    configuration_taus_filename = "results/final_paper_results",simulation_name*"_configuration_taus.csv"
    configuration_betas_filename = "results/final_paper_results",simulation_name*"_configuration_betas.csv"
    energy_taus_filename = "results/final_paper_results",simulation_name*"_energy_taus.csv"
    energy_betas_filename = "results/final_paper_results",simulation_name*"_energy_betas.csv"

    data_matrix = readdlm(joinpath(configuration_taus_filename), ',', Float64, '\n', skipstart=3)
    sample_temperatures = data_matrix[:,1]
    configuration_taus_samples = data_matrix[:,2] 

    data_matrix = readdlm(joinpath(configuration_betas_filename), ',', Float64, '\n', skipstart=3)
    configuration_betas_samples = data_matrix[:,2]

    data_matrix = readdlm(joinpath(energy_taus_filename), ',', Float64, '\n', skipstart=3)
    energy_taus_samples = data_matrix[:,2]

    data_matrix = readdlm(joinpath(energy_betas_filename), ',', Float64, '\n', skipstart=3)
    energy_betas_samples = data_matrix[:,2]

    ### -- CALCULATION ---
    # For all repeated sample tempeartures, remove them from the sample_temperatures array,
    # and collect ALL samples for each temperature from each row and multiple rows and produce a single average per temperature
    # Remember each samples matrix has a row of samples itself
    unique_sample_temperatures = unique(sample_temperatures)
    unique_configuration_taus = zeros(Float64, length(unique_sample_temperatures))
    unique_configuration_betas = zeros(Float64, length(unique_sample_temperatures))
    unique_energy_taus = zeros(Float64, length(unique_sample_temperatures))
    unique_energy_betas = zeros(Float64, length(unique_sample_temperatures))
    for (i, temperature) in pairs(unique_sample_temperatures)
        indices = findall(sample_temperatures .== temperature)
        unique_configuration_taus[i] = mean(configuration_taus_samples[indices])
        unique_configuration_betas[i] = mean(configuration_betas_samples[indices])
        unique_energy_taus[i] = mean(energy_taus_samples[indices])
        unique_energy_betas[i] = mean(energy_betas_samples[indices])
    end

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    Plots.default(dpi = 300)

    ### --- PLOTTING ---
    # Plot all data on same graph, use different y-axis scalings on left and right for taus and BETAS
    graph = plot(title="", xlabel="Temperature, "*L"T", legend=:topleft, yaxis="Relaxation Time, "*L"\tau")
    betas_graph = twinx(graph)
    plot!(graph, unique_sample_temperatures, unique_configuration_taus, label=L"\tau_C", color=alex_red, marker=:circle)
    plot!(graph, unique_sample_temperatures, unique_energy_taus, label=L"\tau_E", color=alex_blue, marker=:circle)

    plot!(betas_graph, unique_sample_temperatures, unique_configuration_betas, label=L"\beta_C", color=alex_orange, marker=:circle,yaxis="Stretching Exponent, "*L"\beta")
    plot!(betas_graph, unique_sample_temperatures, unique_energy_betas, label=L"\beta_E", color=alex_green, marker=:circle)



    ### --- SAVING GRAPHS ---
    savefig(graph, joinpath("results/final_paper_results",simulation_name*"_streched_exponential_parameters_graph.png"))
    savefig(graph, joinpath("results/final_paper_results",simulation_name*"_streched_exponential_parameters_graph.svg"))
    display(graph)




end