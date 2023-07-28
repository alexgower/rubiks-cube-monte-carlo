using DelimitedFiles
using Plots
using StatsBase

function neighbour_energy_deltas_histogram_axis_fixer(data_simulation_name_to_use::String)

    # Read data from file 
    data_matrix = readdlm(joinpath("tools/histogram_axis_fixer/",data_simulation_name_to_use), ',', Any, '\n', skipstart=2)
    neighbour_energy_deltas_sample_temperatures = convert(Vector{Float64},data_matrix[:,1])


    # -------------------------------------

    # New Style
    E_average_by_temperature = convert(Vector{Float64},data_matrix[:,2])

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

    # Find maximum and minima neighbour energy deltas to use as end points of histogram
    max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
    min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)

    # Find maximimum probability density
    max_probability_density = 0
    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
        histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
        probability_densities = collect(values(histogram))
        probability_densities = probability_densities ./ sum(probability_densities)
        this_max_probability_density = max(max_probability_density, maximum(probability_densities))

        if this_max_probability_density > max_probability_density
            max_probability_density = this_max_probability_density
        end
    end

    # Remake histograms ---

    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)

        simulation_name_to_use = "$(data_simulation_name_to_use)_T=$(sample_temperature)"

        # Create a histogram
        histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

        # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
        for energy_delta in min_neighbour_energy_delta:max_neighbour_energy_delta
            if !haskey(histogram, energy_delta)
                histogram[energy_delta] = 0
            end
        end

        energy_deltas = collect(keys(histogram))
        println("Energy deltas: ", energy_deltas)

        probability_densities = collect(values(histogram))
        println("Probability densities: ", probability_densities)

        # Normalize the values
        probability_densities = probability_densities ./ sum(probability_densities)

        # Plot histogram
        graph = bar(energy_deltas, probability_densities, xticks=energy_deltas, xtickfontsize=4, label="Neighbouring Energy Deltas", ylim=(0,max_probability_density*1.1), legend=:topleft, title="Neighbour Energy Deltas Distribution at T=$(sample_temperature)", xlabel="Neighbour Energy Deltas", ylabel="Probability Density")

        # Plot One in a Hundred Energy Dashed Line
        one_in_hundred_energy = log(100)*neighbour_energy_deltas_sample_temperatures[index]
        vline!(graph, [one_in_hundred_energy], linestyle=:dash, label="1/100 Energy", color=:orange)

        # Plot One in a Thousand Energy Dashed Line
        one_in_thousand_energy = log(1000)*neighbour_energy_deltas_sample_temperatures[index]
        vline!(graph, [one_in_thousand_energy], linestyle=:dash, label="1/1000 Energy", color=:red)

        # Save graphs ----------
        savefig(graph, "results/neighbour_energy_deltas_distribution_results/$(simulation_name_to_use).png")

    end

end




function energy_histogram_axis_fixer(data_simulation_name_to_use::String)

    # Read data from file 
    data_matrix = readdlm(joinpath("tools/histogram_axis_fixer/",data_simulation_name_to_use), ',', Any, '\n', skipstart=2)

    energy_sample_temperatures = convert(Vector{Float64},data_matrix[:,1])


    # -------------------------------------


    # # Remove '[' and ']' characters at start and end of neighbour energy delta array
    # start_elements = convert(Vector{String},data_matrix[:,3])
    # start_elements = [element[3:end] for element in start_elements]
    # start_elements = parse.(Float64,start_elements)
    
    # end_elements = convert(Vector{String},data_matrix[:,end])
    # end_elements = [element[1:end-2] for element in end_elements]
    # end_elements = parse.(Float64,end_elements)


    # neighbour_energy_deltas_by_temperature = zeros(length(energy_sample_temperatures), length(data_matrix[1,:])-2) # -2 to remove temperature and E_average columns

    # neighbour_energy_deltas_by_temperature[:,1] .= start_elements
    # neighbour_energy_deltas_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,4:end-1])
    # neighbour_energy_deltas_by_temperature[:,end] .= end_elements

    # OLD STYLE
    # Remove '[' and ']' characters at start and end of neighbour energy delta array
    start_elements = convert(Vector{String},data_matrix[:,2])
    start_elements = [element[3:end] for element in start_elements]
    start_elements = parse.(Float64,start_elements)
    
    end_elements = convert(Vector{String},data_matrix[:,end])
    end_elements = [element[1:end-2] for element in end_elements]
    end_elements = parse.(Float64,end_elements)


    energy_samples_by_temperature = zeros(length(energy_sample_temperatures), length(data_matrix[1,:])-1) # -2 to remove temperature and E_average columns

    energy_samples_by_temperature[:,1] .= start_elements
    energy_samples_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,3:end-1])
    energy_samples_by_temperature[:,end] .= end_elements

    # -------------------------------------

    # Find maximum and minima neighbour energy deltas to use as end points of histogram
    max_energy_sample = maximum(energy_samples_by_temperature)
    min_energy_sample = minimum(energy_samples_by_temperature)

    # Find maximimum probability density
    max_probability_density = 0
    for (index,sample_temperature) in pairs(energy_sample_temperatures)
        histogram = countmap(energy_samples_by_temperature[index,:])
        probability_densities = collect(values(histogram))
        probability_densities = probability_densities ./ sum(probability_densities)
        this_max_probability_density = max(max_probability_density, maximum(probability_densities))

        if this_max_probability_density > max_probability_density
            max_probability_density = this_max_probability_density
        end
    end

    # Remake histograms ---

    for (index,sample_temperature) in pairs(energy_sample_temperatures)

        simulation_name_to_use = "$(data_simulation_name_to_use)_T=$(sample_temperature)"

        # Create a histogram
        histogram = countmap(energy_samples_by_temperature[index,:])

        # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
        for energy_sample_value in min_energy_sample:max_energy_sample
            if !haskey(histogram, energy_sample_value)
                histogram[energy_sample_value] = 0
            end
        end

        energy_sample_values = collect(keys(histogram))
        println("Energy Sample Values: ", energy_sample_values)

        probability_densities = collect(values(histogram))
        println("Probability densities: ", probability_densities)

        # Normalize the values
        probability_densities = probability_densities ./ sum(probability_densities)

        # Plot histogram
        graph = bar(energy_sample_values, probability_densities, xticks=min_energy_sample:20:max_energy_sample, xtickfontsize=4, label="Energy", ylim=(0,max_probability_density), legend=:topleft, xlabel="Energy", ylabel="Probability Density")
        L = simulation_name_to_use[3]
        swap_move_probability = simulation_name_to_use[end-3:end]
        title!(graph, "Energy Histogram for L=$L, P_s=$swap_move_probability, T=$sample_temperature")

        # Save graphs ----------
        savefig(graph, "results/energy_histogram_results/$(simulation_name_to_use).png")

    end

end