using DelimitedFiles
using Plots
using LaTeXStrings
# using StatsBase
# using CSV
# using DataFrames
# using LsqFit
# using Plots.PlotMeasures


include("../core/rubiks_cube.jl")


function phase_transition_graphs_plotter(simulation_stem::String, L_values::Vector{Int64})
    
    temperature_vector_by_L::Vector{Vector{Float64}} = []
    specific_heat_capacities_by_L_by_temperature::Vector{Vector{Float64}} = []
    M_2_by_L_by_temperature::Vector{Vector{Float64}} = []
    M_4_by_L_by_temperature::Vector{Vector{Float64}} = []
    
    
    ##### ----- IMPORT RELAXED ANNEAL DATA -----
    for L in L_values

        filename = joinpath("results/relaxed_anneal_results/"*simulation_stem*"_L="*string(L)*"_1.0")
        # normalization_energy = -solved_configuration_energy(RubiksCube(L))
        normalization_energy = 1.0


        ## -- READ IN THE DATA --
        data_matrix = readdlm(filename, ',', Float64, '\n', skipstart=2, comments=true, comment_char='#')

        # TODO check want normalization energy
        push!(temperature_vector_by_L, copy(data_matrix[:,1]))
        push!(specific_heat_capacities_by_L_by_temperature, copy(data_matrix[:,5]) ./ normalization_energy)
        push!(M_2_by_L_by_temperature, copy(data_matrix[:,6]))
        push!(M_4_by_L_by_temperature, copy(data_matrix[:,7]))


    end


    ### Plot the specific heat capacities against temeprature on the same plot
    specific_heat_capacities_graph = plot(xlabel="Temperature", ylabel="Specific Heat Capacity", title="Specific Heat Capacities with L", legend=:bottomright)
    for L_index in 1:length(L_values)
        plot!(specific_heat_capacities_graph, temperature_vector_by_L[L_index], specific_heat_capacities_by_L_by_temperature[L_index], label="L = "*string(L_values[L_index]))
    end

    savefig(specific_heat_capacities_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_specific_heat_capacities.png"))
    savefig(specific_heat_capacities_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_specific_heat_capacities.svg"))
    display(specific_heat_capacities_graph)

    ###

    ### Make a new figure with maximum specific heat capacities against L

    # Find maximum specific heat capacity for each L_index and their associated temperatures
    max_specific_heat_capacities = []
    max_specific_heat_capacities_temperatures = []
    for L_index in 1:length(L_values)
        max_specific_heat_capacity = maximum(specific_heat_capacities_by_L_by_temperature[L_index])
        max_specific_heat_capacity_temperature = temperature_vector_by_L[L_index][argmax(specific_heat_capacities_by_L_by_temperature[L_index])]
        push!(max_specific_heat_capacities, max_specific_heat_capacity)
        push!(max_specific_heat_capacities_temperatures, max_specific_heat_capacity_temperature)
    end

    max_specific_heat_capacities_graph = plot(xlabel="L", ylabel="Max Specific Heat Capacity", title="Max Specific Heat Capacities with L", legend=false)
    scatter!(max_specific_heat_capacities_graph, L_values, max_specific_heat_capacities, label="Max Specific Heat Capacity")

    savefig(max_specific_heat_capacities_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_max_specific_heat_capacities.png"))
    savefig(max_specific_heat_capacities_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_max_specific_heat_capacities.svg"))
    display(max_specific_heat_capacities_graph)

    ###


    ### Calculate the Binder cumulants for each L against temperature
    binder_cumulants_by_L_by_temperature = []
    for L_index in 1:length(L_values)
        M_2 = M_2_by_L_by_temperature[L_index]
        M_4 = M_4_by_L_by_temperature[L_index]
        push!(binder_cumulants_by_L_by_temperature, 1 .- (M_4 ./ (3 .* M_2.^2)))
    end

    ### Plot the Binder cumulants against temperature on the same plot
    binder_cumulants_graph = plot(xlabel="Temperature", ylabel="Binder Cumulant", title="Binder Cumulants with L", legend=:bottomright)
    for L_index in 1:length(L_values)
        plot!(binder_cumulants_graph, temperature_vector_by_L[L_index], binder_cumulants_by_L_by_temperature[L_index], label="L = "*string(L_values[L_index]))
    end

    savefig(binder_cumulants_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_binder_cumulants.png"))
    savefig(binder_cumulants_graph, joinpath("results/relaxed_anneal_results/"*simulation_stem*"_binder_cumulants.svg"))
    display(binder_cumulants_graph)




end