using DelimitedFiles
using Plots
using StatsBase

function acceptable_neighbours_graph_plotter(graph_save_name::String, data_simulation_name_to_use::String; comparison_data_simulation_name_to_use::String="", log_scale::Bool=false)


    ### ORIGINAL DATA

    # Read data from file 
    data_matrix = readdlm(joinpath("results/tools/acceptable_neighbours_graph_plotter/",data_simulation_name_to_use), ',', Any, '\n', skipstart=2)
    neighbour_energy_deltas_sample_temperatures = convert(Vector{Float64},data_matrix[:,1])


    # -------------------------------------

    # New Style
    # Remove '[' and ']' characters at start and end of neighbour energy delta array
    start_elements = convert(Vector{String},data_matrix[:,3])
    start_elements = [element[3:end] for element in start_elements]
    start_elements = parse.(Float64,start_elements)
    
    end_elements = convert(Vector{String},data_matrix[:,end])
    end_elements = [element[1:end-2] for element in end_elements]
    end_elements = parse.(Float64,end_elements)


    neighbour_energy_deltas_by_temperature = zeros(length(neighbour_energy_deltas_sample_temperatures), length(data_matrix[1,:])-2) # -2 to remove temperature and E_average columns

    neighbour_energy_deltas_by_temperature[:,1] .= start_elements
    neighbour_energy_deltas_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,4:end-1])
    neighbour_energy_deltas_by_temperature[:,end] .= end_elements

    # -------------------------------------


    # OLD STYLE
    # # Remove '[' and ']' characters at start and end of neighbour energy delta array
    # start_elements = convert(Vector{String},data_matrix[:,2])
    # start_elements = [element[3:end] for element in start_elements]
    # start_elements = parse.(Float64,start_elements)
    
    # end_elements = convert(Vector{String},data_matrix[:,end])
    # end_elements = [element[1:end-2] for element in end_elements]
    # end_elements = parse.(Float64,end_elements)
    

    # neighbour_energy_deltas_by_temperature = zeros(length(neighbour_energy_deltas_sample_temperatures), length(data_matrix[1,:])-1) # -2 to remove temperature and E_average columns

    # neighbour_energy_deltas_by_temperature[:,1] .= start_elements
    # neighbour_energy_deltas_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,3:end-1])
    # neighbour_energy_deltas_by_temperature[:,end] .= end_elements

    # -------------------------------------



    ### COMPARISON DATA

    if comparison_data_simulation_name_to_use != ""

        # Read data from file 
        data_matrix = readdlm(joinpath("results/tools/acceptable_neighbours_graph_plotter/",comparison_data_simulation_name_to_use), ',', Any, '\n', skipstart=2)
        comparison_neighbour_energy_deltas_sample_temperatures = convert(Vector{Float64},data_matrix[:,1])


        # -------------------------------------

        # New Style
        # Remove '[' and ']' characters at start and end of neighbour energy delta array
        start_elements = convert(Vector{String},data_matrix[:,3])
        start_elements = [element[3:end] for element in start_elements]
        start_elements = parse.(Float64,start_elements)
        
        end_elements = convert(Vector{String},data_matrix[:,end])
        end_elements = [element[1:end-2] for element in end_elements]
        end_elements = parse.(Float64,end_elements)


        comparison_neighbour_energy_deltas_by_temperature = zeros(length(comparison_neighbour_energy_deltas_sample_temperatures), length(data_matrix[1,:])-2) # -2 to remove temperature and E_average columns

        comparison_neighbour_energy_deltas_by_temperature[:,1] .= start_elements
        comparison_neighbour_energy_deltas_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,4:end-1])
        comparison_neighbour_energy_deltas_by_temperature[:,end] .= end_elements
        
        # -------------------------------------

    end


    ### ORIGINAL DATA

    # Calculate Acceptable Neighbour Proportions

    one_in_hundred_acceptable_neighbour_proportions_by_temperature = zeros(length(neighbour_energy_deltas_sample_temperatures))
    one_in_thousand_acceptable_neighbour_proportions_by_temperature = zeros(length(neighbour_energy_deltas_sample_temperatures))


    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)

        # Create a histogram
        histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

        energy_deltas = collect(keys(histogram))
        println("Energy deltas: ", energy_deltas)

        probability_densities = collect(values(histogram))
        println("Probability densities: ", probability_densities)

        # Normalize the values
        probability_densities = probability_densities ./ sum(probability_densities)

        # Calculate probability density weight below One in a Hundred Energy Dashed Line
        one_in_hundred_energy = log(100)*neighbour_energy_deltas_sample_temperatures[index]
        one_in_hundred_energy_index = findfirst(energy_deltas .> one_in_hundred_energy)

        if one_in_hundred_energy_index == nothing
            one_in_hundred_energy_index = length(energy_deltas)
        end

        probability_density_weight_below_one_in_hundred_energy = sum(probability_densities[1:one_in_hundred_energy_index])
        one_in_hundred_acceptable_neighbour_proportions_by_temperature[index] = probability_density_weight_below_one_in_hundred_energy

        # Calculate probability density weight below One in a Thousand Energy Dashed Line
        one_in_thousand_energy = log(1000)*neighbour_energy_deltas_sample_temperatures[index]
        one_in_thousand_energy_index = findfirst(energy_deltas .> one_in_thousand_energy)

        if one_in_thousand_energy_index == nothing
            one_in_thousand_energy_index = length(energy_deltas)
        end
        
        probability_density_weight_below_one_in_thousand_energy = sum(probability_densities[1:one_in_thousand_energy_index])
        one_in_thousand_acceptable_neighbour_proportions_by_temperature[index] = probability_density_weight_below_one_in_thousand_energy

    end


    graph = log_scale ? plot(xlabel="Temperature", ylabel="Acceptable Neighbour Proportion", title="Acceptable Neighbour Proportion vs Temperature", xscale=:log10, yscale=:log10, legend=:topright) : plot(xlabel="Temperature", ylabel="Acceptable Neighbour Proportion", title="Acceptable Neighbour Proportion vs Temperature", legend=:bottomleft, ylimits=(0,1))
   
    ### ADD T_c and T^* LINES TO GRAPH
    T_c = 0.955
    T_star = 1.0
    T_star_low = 0.93
    T_star_end = 0.84
    vline!(graph, [T_c], linestyle=:dash, label="T_c", color=:black, lineopacity=0.5)
    vline!(graph, [T_star], linestyle=:dash, label="T^*", color=:purple, lineopacity=0.5)
    # vline!(graph, [T_star_low], linestyle=:dash, label="T^* Low", color=:orange)
    # vline!(graph, [T_star_end], linestyle=:dash, label="T^* End", color=:purple)

   
   
   
    # Plot Data
    plot!(graph, neighbour_energy_deltas_sample_temperatures, one_in_hundred_acceptable_neighbour_proportions_by_temperature, label="1/100 Energy (6 Swap Moves)", color=:orange)
    plot!(graph, neighbour_energy_deltas_sample_temperatures, one_in_thousand_acceptable_neighbour_proportions_by_temperature, label="1/1000 Energy (6 Swap Move)", color=:red)




    ### COMPARISON

    if comparison_data_simulation_name_to_use != ""

        # Calculate Acceptable Neighbour Proportions

        comparison_one_in_hundred_acceptable_neighbour_proportions_by_temperature = zeros(length(comparison_neighbour_energy_deltas_sample_temperatures))
        comparison_one_in_thousand_acceptable_neighbour_proportions_by_temperature = zeros(length(comparison_neighbour_energy_deltas_sample_temperatures))


        for (index,sample_temperature) in pairs(comparison_neighbour_energy_deltas_sample_temperatures)

            # Create a histogram
            histogram = countmap(comparison_neighbour_energy_deltas_by_temperature[index,:])

            energy_deltas = collect(keys(histogram))
            println("Energy deltas: ", energy_deltas)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Calculate probability density weight below One in a Hundred Energy Dashed Line
            one_in_hundred_energy = log(100)*comparison_neighbour_energy_deltas_sample_temperatures[index]
            one_in_hundred_energy_index = findfirst(energy_deltas .> one_in_hundred_energy)

            if one_in_hundred_energy_index == nothing
                one_in_hundred_energy_index = length(energy_deltas)
            end

            probability_density_weight_below_one_in_hundred_energy = sum(probability_densities[1:one_in_hundred_energy_index])
            comparison_one_in_hundred_acceptable_neighbour_proportions_by_temperature[index] = probability_density_weight_below_one_in_hundred_energy

            # Calculate probability density weight below One in a Thousand Energy Dashed Line
            one_in_thousand_energy = log(1000)*comparison_neighbour_energy_deltas_sample_temperatures[index]
            one_in_thousand_energy_index = findfirst(energy_deltas .> one_in_thousand_energy)

            if one_in_thousand_energy_index == nothing
                one_in_thousand_energy_index = length(energy_deltas)
            end
            
            probability_density_weight_below_one_in_thousand_energy = sum(probability_densities[1:one_in_thousand_energy_index])
            comparison_one_in_thousand_acceptable_neighbour_proportions_by_temperature[index] = probability_density_weight_below_one_in_thousand_energy

        end


        # Plot Graph
        plot!(graph, comparison_neighbour_energy_deltas_sample_temperatures, comparison_one_in_hundred_acceptable_neighbour_proportions_by_temperature, label="1/100 Energy (1 Swap Move)", color=:green)
        plot!(graph, comparison_neighbour_energy_deltas_sample_temperatures, comparison_one_in_thousand_acceptable_neighbour_proportions_by_temperature, label="1/1000 Energy (1 Swap Move)", color=:blue)
            
    end




    
    ### FLIP GRAPH
    xflip!(graph)

    ### ADJUST LEGEND
    graph = plot!(legend=:bottomleft, legendfontsize=15)

    ### SAVE GRAPH
    savefig(graph, "results/tools/acceptable_neighbours_graph_plotter/$(graph_save_name).png")


    printstyled("Original Data 1/1000 Proportions: $(reverse(one_in_thousand_acceptable_neighbour_proportions_by_temperature))"; color=:light_blue)
    if comparison_data_simulation_name_to_use != ""
        printstyled("Comparison Data 1/1000 Proportions: $(reverse(comparison_one_in_thousand_acceptable_neighbour_proportions_by_temperature))"; color=:light_red)
    end

    display(graph)
end
